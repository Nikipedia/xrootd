//------------------------------------------------------------------------------
// Copyright (c) 2011-2014 by European Organization for Nuclear Research (CERN)
// Author: Michal Simon <michal.simon@cern.ch>
//------------------------------------------------------------------------------
// This file is part of the XRootD software suite.
//
// XRootD is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// XRootD is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with XRootD.  If not, see <http://www.gnu.org/licenses/>.
//
// In applying this licence, CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
//------------------------------------------------------------------------------

#include "XrdEc/XrdEcReader.hh"
#include "XrdEc/XrdEcUtilities.hh"
#include "XrdEc/XrdEcConfig.hh"
#include "XrdEc/XrdEcObjCfg.hh"
#include "XrdEc/XrdEcThreadPool.hh"

#include "XrdZip/XrdZipCDFH.hh"
#include "XrdZip/XrdZipLFH.hh"
#include "XrdZip/XrdZipUtils.hh"

#include "XrdCl/XrdClMessageUtils.hh"


#include "XrdOuc/XrdOucCRC32C.hh"

#include "XrdCl/XrdClParallelOperation.hh"
#include "XrdCl/XrdClZipOperations.hh"
#include "XrdCl/XrdClFileOperations.hh"
#include "XrdCl/XrdClFinalOperation.hh"

#include <algorithm>
#include <iterator>
#include <numeric>
#include <tuple>
#include <iostream>

namespace XrdEc
{
  //---------------------------------------------------------------------------
  // Destructor (we need it in the source file because block_t is defined in
  // here)
  //---------------------------------------------------------------------------
  Reader::~Reader()
  {
  }

  //---------------------------------------------------------------------------
  // Open the erasure coded / striped object
  //---------------------------------------------------------------------------
  void Reader::Open( XrdCl::ResponseHandler *handler, uint16_t timeout )
  {
    const size_t size = objcfg.plgr.size();
    std::vector<XrdCl::Pipeline> opens; opens.reserve( size );
    std::vector<XrdCl::Pipeline> healthRead; healthRead.reserve(size);
    for( size_t i = 0; i < size; ++i )
    {
      // generate the URL
      std::string url = objcfg.GetDataUrl( i );

      // create the file object
      dataarchs.emplace( url, std::make_shared<XrdCl::ZipArchive>(
          Config::Instance().enable_plugins ) );
      // open the archive
      if( objcfg.nomtfile )
      {
        opens.emplace_back( XrdCl::OpenArchive( *dataarchs[url], url, XrdCl::OpenFlags::Read ) );
      }
      else
        opens.emplace_back( OpenOnly( *dataarchs[url], url, false ) );

      healthRead.emplace_back(ReadHealth(i));
    }



    auto pipehndl = [=]( const XrdCl::XRootDStatus &st )
                    { // set the central directories in ZIP archives (if we use metadata files)
        auto itr = dataarchs.begin();
		for (; itr != dataarchs.end(); ++itr) {
			const std::string &url = itr->first;
			auto &zipptr = itr->second;
			if (zipptr->openstage == XrdCl::ZipArchive::NotParsed)
				zipptr->SetCD(metadata[url]);
			else if (zipptr->openstage != XrdCl::ZipArchive::Done
					&& !metadata.empty())
				AddMissing(metadata[url]);
			auto itr = zipptr->cdmap.begin();
			for (; itr != zipptr->cdmap.end(); ++itr) {
				try {
					size_t blknb = fntoblk(itr->first);
					urlmap.emplace(itr->first, url);
					if (blknb > lstblk)
						lstblk = blknb;
				} catch (std::invalid_argument&) {
					std::cout << "Invalid file name detected\n" << std::flush;
				}
			}
		}
                      metadata.clear();
                      // call user handler
                      if( handler )
                        handler->HandleResponse( new XrdCl::XRootDStatus( st ), nullptr );
                    };
    // in parallel open the data files and read the metadata
    XrdCl::Pipeline p = objcfg.nomtfile
                      ? XrdCl::Parallel( opens ).AtLeast( objcfg.nbdata )
                    		  //| XrdCl::Parallel(healthRead).AtLeast(objcfg.nbdata)
                    		  | ReadSize( 0 ) | XrdCl::Final( pipehndl )
                      : XrdCl::Parallel( ReadMetadata( 0 ),
                                         XrdCl::Parallel( opens ).AtLeast( objcfg.nbdata ) ) >> pipehndl;
    XrdCl::Async( std::move( p ), timeout );
  }

  //-----------------------------------------------------------------------
  // Read data from the data object
  //-----------------------------------------------------------------------
  void Reader::Read( uint64_t                offset,
                     uint32_t                length,
                     void                   *buffer,
                     XrdCl::ResponseHandler *handler,
                     uint16_t                timeout,
					 bool 					 doRepair)
  {
    if( objcfg.nomtfile )
    {
    	if(offset >= filesize){
    		length = 0;
    	}
    	else if( offset + length > filesize )
        {
    	  length = filesize - offset;
        }
    }

    if(length == 0){
    	ScheduleHandler( offset, 0, buffer, handler);
    	return;
    }

    char *usrbuff = reinterpret_cast<char*>( buffer );
    typedef std::tuple<uint64_t, uint32_t,
                       void*, uint32_t,
                       XrdCl::ResponseHandler*,
                       XrdCl::XRootDStatus> rdctx_t;
    auto rdctx = std::make_shared<rdctx_t>( offset, 0, buffer,
                                            length, handler,
                                            XrdCl::XRootDStatus() );
    auto rdmtx = std::make_shared<std::mutex>();

    while( length > 0 )
    {
      size_t   blkid  = offset / objcfg.datasize;                                     //< ID of the block from which we will be reading
      size_t   strpid = ( offset % objcfg.datasize ) / objcfg.chunksize;              //< ID of the stripe from which we will be reading
      uint64_t rdoff  = offset - blkid * objcfg.datasize - strpid * objcfg.chunksize; //< relative read offset within the stripe
      uint32_t rdsize = objcfg.chunksize - rdoff;                                     //< read size within the stripe
      if( rdsize > length ) rdsize = length;
      //-------------------------------------------------------------------
      // Make sure we operate on a valid block
      //-------------------------------------------------------------------
      std::unique_lock<std::mutex> lck( blkmtx );
      if( !block || block->blkid != blkid )
        block = std::make_shared<block_t>( blkid, *this, objcfg );
      //-------------------------------------------------------------------
      // Prepare the callback for reading from single stripe
      //-------------------------------------------------------------------
      auto blk = block;
      lck.unlock();
      auto callback = [blk, rdctx, rdsize, rdmtx]( const XrdCl::XRootDStatus &st, uint32_t nbrd )
      {
    	  std::unique_lock<std::mutex> lck( *rdmtx );
        //---------------------------------------------------------------------
        // update number of bytes left to be read (bytes requested not actually
        // read)
        //---------------------------------------------------------------------
        std::get<3>( *rdctx ) -= rdsize;
        //---------------------------------------------------------------------
        // Handle failure ...
        //---------------------------------------------------------------------
        if( !st.IsOK() )
          std::get<5>( *rdctx ) = st; // the error
        //---------------------------------------------------------------------
        // Handle success ...
        //---------------------------------------------------------------------
        else
          std::get<1>( *rdctx ) += nbrd; // number of bytes read
        //---------------------------------------------------------------------
        // Are we done?
        //---------------------------------------------------------------------
        if( std::get<3>( *rdctx ) == 0 )
        {
          //-------------------------------------------------------------------
          // Check if the read operation was successful ...
          //-------------------------------------------------------------------
          XrdCl::XRootDStatus &status = std::get<5>( *rdctx );
          if( !status.IsOK() )
            ScheduleHandler( std::get<4>( *rdctx ), status );
          else
            ScheduleHandler( std::get<0>( *rdctx ), std::get<1>( *rdctx ),
                             std::get<2>( *rdctx ), std::get<4>( *rdctx ) );
        }
      };
      //-------------------------------------------------------------------
      // Read data from a stripe
      //-------------------------------------------------------------------
      block_t::read( blk, strpid, rdoff, rdsize, usrbuff, callback, timeout, doRepair );
      //-------------------------------------------------------------------
      // Update absolute offset, read length, and user buffer
      //-------------------------------------------------------------------
      offset  += rdsize;
      length  -= rdsize;
      usrbuff += rdsize;
    }
  }

  //-----------------------------------------------------------------------
  // Close the data object
  //-----------------------------------------------------------------------
  void Reader::Close( XrdCl::ResponseHandler *handler, uint16_t timeout )
  {
    //---------------------------------------------------------------------
    // prepare the pipelines ...
    //---------------------------------------------------------------------
    std::vector<XrdCl::Pipeline> closes;
    closes.reserve( dataarchs.size() );
    auto itr = dataarchs.begin();
    for( ; itr != dataarchs.end() ; ++itr )
    {
      auto &zipptr = itr->second;
      if( zipptr->IsOpen() )
      {
        zipptr->SetProperty( "BundledClose", "true");
        closes.emplace_back( XrdCl::CloseArchive( *zipptr ) >>
                             [zipptr]( XrdCl::XRootDStatus& ){ } );
      }
    }

    // if there is nothing to close just schedule the handler
    if( closes.empty() ) ScheduleHandler( handler );
    // otherwise close the archives
    else XrdCl::Async( XrdCl::Parallel( closes ) >> handler, timeout );
  }

  void Reader::DebuggingRead(size_t blknb, size_t strpnb, buffer_t &buffer, callback_t cb, uint16_t timeout)
  {
	  Read(blknb, strpnb, buffer, cb, timeout, true);
  }

  //-------------------------------------------------------------------------
  // on-definition is not allowed here beforeiven stripes from given block
  //-------------------------------------------------------------------------
  void Reader::Read( size_t blknb, size_t strpnb, buffer_t &buffer, callback_t cb, uint16_t timeout, bool exactControl )
  {
    // generate the file name (blknb/strpnb)
    std::string fn = objcfg.GetFileName( blknb, strpnb );
    // if the block/stripe does not exist it means we are reading passed the end of the file
    auto itr = urlmap.find( fn );
    if( itr == urlmap.end() )
    {
      //auto st = !IsMissing( fn ) ? XrdCl::XRootDStatus() :
      //          XrdCl::XRootDStatus( XrdCl::stError, XrdCl::errNotFound );
    	auto st = XrdCl::XRootDStatus( XrdCl::stError, XrdCl::errNotFound );
    	ThreadPool::Instance().Execute( cb, st, 0 );
      return;
    }
    // get the URL of the ZIP archive with the respective data
    const std::string &url = itr->second;
    // get the ZipArchive object
    auto &zipptr = dataarchs[url];
    // check the size of the data to be read
    XrdCl::StatInfo *info = nullptr;
    auto st = zipptr->Stat( fn, info );
    if( !st.IsOK() )
    {
    	ThreadPool::Instance().Execute( cb, st, 0 );
      return;
    }
    uint32_t rdsize = info->GetSize();
    delete info;
    // create a buffer for the data
    buffer.resize( objcfg.chunksize );
    // issue the read request
    XrdCl::Async(XrdCl::ReadFrom( *zipptr, fn, 0, rdsize, buffer.data() ) >>
                    [zipptr, fn, cb, &buffer, exactControl, this, url, timeout]( XrdCl::XRootDStatus &st, XrdCl::ChunkInfo &ch )
                    {
                      //---------------------------------------------------
                      // If read failed there's nothing to do, just pass the
                      // status to user callback
                      //---------------------------------------------------
                      if( !st.IsOK() )
                      {
                    	  cb( XrdCl::XRootDStatus(st.status, "Read failed"), 0 );
                        return;
                      }
                      //---------------------------------------------------
                      // Get the checksum for the read data
                      //---------------------------------------------------
                      uint32_t orgcksum = 0;
                      //auto s = zipptr->GetCRC32( fn, orgcksum );
                      auto s = zipptr->GetCRC32(fn, orgcksum);
                      //---------------------------------------------------
                      // If we cannot extract the checksum assume the data
                      // are corrupted
                      //---------------------------------------------------
                      if( !s.IsOK() )
                      {
                        cb( XrdCl::XRootDStatus(s.status, s.code, s.errNo, "Chksum fail"), 0 );
                        return;
                      }
                      //---------------------------------------------------
                      // Verify data integrity
                      //---------------------------------------------------
                      uint32_t cksum = objcfg.digest( 0, ch.buffer, ch.length );
                      if( orgcksum != cksum )
                      {
                    	  cb( XrdCl::XRootDStatus( XrdCl::stError, "Chksum unequal" ), 0 );
                        return;
                      }
                      // optionally also checks LFH against CDFH
                      if (exactControl) {
                    	  auto cditr = zipptr->cdmap.find(fn);
							if (cditr == zipptr->cdmap.end())
								{
								cb(XrdCl::XRootDStatus(XrdCl::stError, "File not found."),0);
								return;
								}


							XrdZip::CDFH *cdfh =
									zipptr->cdvec[cditr->second].get();
							uint32_t offset = XrdZip::CDFH::GetOffset(*cdfh);
							uint64_t cdOffset =
									zipptr->zip64eocd ?
											zipptr->zip64eocd->cdOffset :
											zipptr->eocd->cdOffset;
							uint64_t nextRecordOffset =
									(cditr->second + 1 < zipptr->cdvec.size()) ?
											XrdZip::CDFH::GetOffset(
													*zipptr->cdvec[cditr->second
															+ 1]) :
											cdOffset;
							auto readSize = (nextRecordOffset - offset)- cdfh->uncompressedSize;
							std::shared_ptr<buffer_t> lfhbuf;
							lfhbuf = std::make_shared<buffer_t>();
							lfhbuf->reserve(readSize);
							XrdCl::Pipeline p =  XrdCl::Read(zipptr->archive, offset, readSize, lfhbuf->data(), timeout) >>
							                   [=]( XrdCl::XRootDStatus &st, XrdCl::ChunkInfo &ch )
							                   {
												if (st.IsOK()) {
													try {
														XrdZip::LFH lfh(
																lfhbuf->data());
												if (lfh.ZCRC32 != orgcksum) {
													cb(
															XrdCl::XRootDStatus(
																	XrdCl::stError,
	"CRC!=CDFH"),
	0);
	return;
}
												if (lfh.compressedSize
														!= cdfh->compressedSize) {
													std::stringstream ss;
													ss << "CompSize " << lfh.compressedSize << " != " << cdfh->compressedSize;
													cb(
															XrdCl::XRootDStatus(
																	XrdCl::stError,
																	ss.str()),
															0);
													return;
												}
												if (lfh.compressionMethod
														!= cdfh->compressionMethod
														|| lfh.extraLength
																!= cdfh->extraLength){
													std::stringstream ss;
													ss << "CompMethod "
															<< lfh.compressionMethod
															<< " != "
															<< cdfh->compressionMethod
															<< " extraLength "
															<< lfh.extraLength
															<< " != "
															<< cdfh->extraLength;
													cb(
															XrdCl::XRootDStatus(
																	XrdCl::stError,
																	ss.str()),
															0);
													return;
												}
														if( lfh.filename
																!= cdfh->filename
														|| lfh.filenameLength
																!= cdfh->filenameLength){
															std::stringstream ss;
															ss << "filename!=CDFH" << lfh.filename << ":"<<cdfh->filename<< " lengths: " << lfh.filenameLength << ":" << cdfh->filenameLength << "\n" << std::flush;
															cb(
																	XrdCl::XRootDStatus(
																			XrdCl::stError,
																			ss.str()),
																	0);
															return;
														}
														if( lfh.generalBitFlag
																!= cdfh->generalBitFlag
														|| lfh.minZipVersion
																!= cdfh->minZipVersion) {
													cb(
															XrdCl::XRootDStatus(
																	XrdCl::stError,
																	"LFH!=CDFH"),
															0);
													return;
												}
												if (lfh.uncompressedSize
														!= cdfh->uncompressedSize) {
													cb(
															XrdCl::XRootDStatus(
																	XrdCl::stError,
																	"uncompSize!=CDFH"),
															0);
													return;
												}
														//---------------------------------------------------
														// All is good, we can call now the user callback
														//---------------------------------------------------
														else {
															cb(
																	XrdCl::XRootDStatus(),
																	ch.length);
															return;
														}
													} catch (const std::exception &e) {

													}
													cb(
															XrdCl::XRootDStatus(
																	XrdCl::stError,
																	XrdCl::errCorruptedHeader,
																	0,
																	"Couldnt parse lfh"),
															0);
													return;
												} else {
													cb(XrdCl::XRootDStatus(st.status, st.code, st.errNo, "ReadFrom err"), 0);
																					return;
												}
											};
							    Async( std::move( p ), timeout );
							    return;
						} else {
							//---------------------------------------------------
							// All is good, we can call now the user callback
							//---------------------------------------------------
							cb(XrdCl::XRootDStatus(), ch.length);
							return;
						}
                    }, timeout);
    }

  //-----------------------------------------------------------------------
  // Read metadata for the object
  //-----------------------------------------------------------------------
  XrdCl::Pipeline Reader::ReadMetadata( size_t index )
  {
    const size_t size = objcfg.plgr.size();
    // create the File object
    auto file = std::make_shared<XrdCl::File>( Config::Instance().enable_plugins );
    // prepare the URL for Open operation
    std::string url = objcfg.GetMetadataUrl( index );
    // arguments for the Read operation
    XrdCl::Fwd<uint32_t> rdsize;
    XrdCl::Fwd<void*>    rdbuff;

    return XrdCl::Open( *file, url, XrdCl::OpenFlags::Read ) >>
             [=]( XrdCl::XRootDStatus &st, XrdCl::StatInfo &info ) mutable
             {
               if( !st.IsOK() )
               {
                 if( index + 1 < size )
                   XrdCl::Pipeline::Replace( ReadMetadata( index + 1 ) );
                 return;
               }
               // prepare the args for the subsequent operation
               rdsize = info.GetSize();
               rdbuff = new char[info.GetSize()];
             }
         | XrdCl::Read( *file, 0, rdsize, rdbuff ) >>
             [=]( XrdCl::XRootDStatus &st, XrdCl::ChunkInfo &ch )
             {
               if( !st.IsOK() )
               {
                 if( index + 1 < size )
                   XrdCl::Pipeline::Replace( ReadMetadata( index + 1 ) );
                 return;
               }
               // now parse the metadata
               if( !ParseMetadata( ch ) )
               {
                 if( index + 1 < size )
                   XrdCl::Pipeline::Replace( ReadMetadata( index + 1 ) );
                 return;
               }
             }
         | XrdCl::Close( *file ) >>
             []( XrdCl::XRootDStatus &st )
             {
               if( !st.IsOK() )
                 XrdCl::Pipeline::Ignore(); // ignore errors, we don't really care
             }
         | XrdCl::Final(
             [rdbuff, file]( const XrdCl::XRootDStatus& )
             {
               // deallocate the buffer if necessary
               if( rdbuff.Valid() )
               {
                 char* buffer = reinterpret_cast<char*>( *rdbuff );
                 delete[] buffer;
               }
             } );
  }

  //-----------------------------------------------------------------------
  //! Read size from xattr
  //!
  //! @param index : placement's index
  //-----------------------------------------------------------------------
  XrdCl::Pipeline Reader::ReadSize( size_t index )
  {
    std::string url = objcfg.GetDataUrl( index );
    return XrdCl::GetXAttr( dataarchs[url]->GetFile(), "xrdec.filesize" ) >>
        [index, this, url]( XrdCl::XRootDStatus &st, std::string &size)
        {
          if( !st.IsOK() || this->dataarchs[url]->openstage != XrdCl::ZipArchive::Done)
          {
            //-------------------------------------------------------------
            // Check if we can recover the error or a diffrent location
            //-------------------------------------------------------------
            if( index + 1 < objcfg.plgr.size() )
              XrdCl::Pipeline::Replace( ReadSize( index + 1 ) );
            return;
          }
          try{
          filesize = std::stoull( size );
          }
          catch(std::invalid_argument &){
        	  if( index + 1 < objcfg.plgr.size() )
        	  XrdCl::Pipeline::Replace( ReadSize( index + 1 ) );
          }
        };
  }

  XrdCl::Pipeline Reader::ReadHealth(size_t index){
		  std::string url = objcfg.GetDataUrl( index );
		  return XrdCl::GetXAttr( dataarchs[url]->GetFile(), "xrdec.unhealthy" ) >>
	          [index, url, this]( XrdCl::XRootDStatus &st, std::string &damage)
	          {
	            if( !st.IsOK() )
	            {
	              //-------------------------------------------------------------
	              // Check if we can recover the error or a diffrent location
	              //-------------------------------------------------------------

	              return;
	            }
	            int damaged = std::stoi( damage );
	            if(damaged > 0) this->dataarchs[url]-> openstage = XrdCl::ZipArchive::Error;
	          };
  }

  //-----------------------------------------------------------------------
  // Parse metadata from chunk info object
  //-----------------------------------------------------------------------
  bool Reader::ParseMetadata( XrdCl::ChunkInfo &ch )
  {
    const size_t mincnt = objcfg.nbdata + objcfg.nbparity;
    const size_t maxcnt = objcfg.plgr.size();

    char   *buffer = reinterpret_cast<char*>( ch.buffer );
    size_t  length = ch.length;

    for( size_t i = 0; i < maxcnt; ++i )
    {
      uint32_t signature = XrdZip::to<uint32_t>( buffer );
      if( signature != XrdZip::LFH::lfhSign )
      {
        if( i + 1 < mincnt ) return false;
        break;
      }
      XrdZip::LFH lfh( buffer );
      // check if we are not reading passed the end of the buffer
      if( lfh.lfhSize + lfh.uncompressedSize > length ) return false;
      buffer += lfh.lfhSize;
      length -= lfh.lfhSize;
      // verify the checksum
      uint32_t crc32val = objcfg.digest( 0, buffer, lfh.uncompressedSize );
      if( crc32val != lfh.ZCRC32 ) return false;
      // keep the metadata
      try{
    	  std::string url = objcfg.GetDataUrl( std::stoull( lfh.filename ) );
    	  metadata.emplace( url, buffer_t( buffer, buffer + lfh.uncompressedSize ) );
      }
      catch(const std::invalid_argument* e){
    	  std::cout << "Invalid filename found";
      }
      buffer += lfh.uncompressedSize;
      length -= lfh.uncompressedSize;
    }

    return true;
  }

  //-----------------------------------------------------------------------
  // Add all the entries from given Central Directory to missing
  //-----------------------------------------------------------------------
  void Reader::AddMissing( const buffer_t &cdbuff )
  {
    const char *buff = cdbuff.data();
    size_t      size = cdbuff.size();
    // parse Central Directory records
    XrdZip::cdvec_t cdvec;
    XrdZip::cdmap_t cdmap;
    std::tie(cdvec, cdmap ) = XrdZip::CDFH::Parse( buff, size );
    auto itr = cdvec.begin();
    for( ; itr != cdvec.end() ; ++itr )
    {
      XrdZip::CDFH &cdfh = **itr;
      missing.insert( cdfh.filename );
    }
  }

  //-----------------------------------------------------------------------
  //! Check if chunk file name is missing
  //-----------------------------------------------------------------------
  bool Reader::IsMissing( const std::string &fn )
  {
    // if the chunk is in the missing set return true
    if( missing.count( fn ) ) return true;
    // if we don't have a metadata file and the chunk exceeds last chunk
    // also return true
    try{
    if( objcfg.nomtfile && fntoblk( fn ) <= lstblk ) return true;}
    catch(...){}
    // otherwise return false
    return false;
  }

} /* namespace XrdEc */
