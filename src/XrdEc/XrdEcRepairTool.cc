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

#include "XrdEcRepairTool.hh"

#include "XrdEc/XrdEcReader.hh"
#include "XrdEc/XrdEcUtilities.hh"
#include "XrdEc/XrdEcConfig.hh"
#include "XrdEc/XrdEcObjCfg.hh"
#include "XrdEc/XrdEcThreadPool.hh"

#include "XrdZip/XrdZipLFH.hh"
#include "XrdZip/XrdZipCDFH.hh"
#include "XrdZip/XrdZipEOCD.hh"
#include "XrdZip/XrdZipUtils.hh"

#include "XrdOuc/XrdOucCRC32C.hh"

#include "XrdCl/XrdClParallelOperation.hh"
#include "XrdCl/XrdClZipOperations.hh"
#include "XrdCl/XrdClFileOperations.hh"
#include "XrdCl/XrdClFinalOperation.hh"
#include "XrdCl/XrdClMessageUtils.hh"

#include <algorithm>
#include <iterator>
#include <numeric>
#include <tuple>
#include <iostream>
#include <fstream>

namespace XrdEc {


//-----------------------------------------------------------------------
// If neccessary trigger error correction procedure
// @param self : the block_t object
// @return     : false if the block is corrupted and cannot be recovered,
//               true otherwise
//-----------------------------------------------------------------------
bool RepairTool::error_correction( std::shared_ptr<block_t> &self, RepairTool *writer )
{

    //---------------------------------------------------------------------
    // Do the accounting for our stripes
    //---------------------------------------------------------------------
    size_t missingcnt = 0, validcnt = 0, loadingcnt = 0, recoveringcnt = 0, emptycnt = 0;
    std::for_each( self->state.begin(), self->state.end(), [&]( block_t::state_t &s )
      {
        switch( s )
        {
          case block_t::Missing:    ++missingcnt;    break;
          case block_t::Valid:      ++validcnt;      break;
          case block_t::Loading:    ++loadingcnt;    break;
          case block_t::Recovering: ++recoveringcnt; break;
          case block_t::Empty: 		++emptycnt; 	 break;
          default: ;
        }
      } );

	if(validcnt == writer->objcfg.nbchunks){
    	// both check_block and update_callback will skip to the next block by returning false.
    	return false;
    }
    //---------------------------------------------------------------------
    // Check if we can do the recovery at all (if too many stripes are
    // missing it won't be possible)
    //---------------------------------------------------------------------
    if( missingcnt + recoveringcnt > self->objcfg.nbparity )
    {
    	std::cout<<"Recovery impossible\n"<<std::flush;
      std::for_each( self->state.begin(), self->state.end(),
                     []( block_t::state_t &s ){ if( s == block_t::Recovering ) s = block_t::Missing; } );
      writer->repairFailed = true;
      return false;
    }
    //---------------------------------------------------------------------
    // Check if we can do the recovery right away
    //---------------------------------------------------------------------
    if( validcnt >= self->objcfg.nbdata && missingcnt > 0 )
    {

      Config &cfg = Config::Instance();
      stripes_t strps( self->get_stripes() );
      try
      {
        cfg.GetRedundancy( self->objcfg ).compute( strps );
      }
      catch( const IOError &ex )
      {
        std::for_each( self->state.begin(), self->state.end(),
                       []( block_t::state_t &s ){ if( s == block_t::Recovering ) s = block_t::Missing; } );
        writer->repairFailed = true;
        return false;
      }
      //-------------------------------------------------------------------
      // Now that we triggered the recovery procedure mark every missing
      // stripe as recovering.
      //-------------------------------------------------------------------
      std::for_each( self->state.begin(), self->state.end(),
                     []( block_t::state_t &s ){ if( s == block_t::Missing ) s = block_t::Recovering; } );

      //-------------------------------------------------------------------
      // Now when we recovered the data we need to mark every stripe as
      // valid and execute the pending reads
      //-------------------------------------------------------------------
      for( size_t strpid = 0; strpid < self->objcfg.nbchunks; ++strpid )
      {
        if( self->state[strpid] != block_t::Recovering ) continue;
        // Write new content to disk
        auto st = writer->WriteChunk(self, strpid);
        if(!st.IsOK()){
        	std::cout<<"Writing to archive failed: " << st.GetErrorMessage()<< "\n"<<std::flush;
        	      std::for_each( self->state.begin(), self->state.end(),
        	                     []( block_t::state_t &s ){ if( s == block_t::Recovering ) s = block_t::Missing; } );
        	      writer->repairFailed = true;
        	      return false;
        }
        // TODO: Check this code
        if(writer->checkAfterRepair){
        	writer->Read( self->blkid, strpid, self->stripes[strpid],
        	                   RepairTool::update_callback( self, writer, strpid ), 0, true );
        }
        else{
        	writer->chunksRepaired.fetch_add(1);
        	self->state[strpid] = block_t::Valid;
        	validcnt++;
        	if(validcnt == writer->objcfg.nbchunks){
        			std::cout<<"All Errors corrected\n"<<std::flush;
        	    	return false;
        	    }
        }
      }
      return true;
    }
    //---------------------------------------------------------------------
    // Try loading the data and only then attempt recovery
    //---------------------------------------------------------------------
    size_t i = 0;
    while( loadingcnt + validcnt < self->objcfg.nbchunks && i < self->objcfg.nbchunks )
    {
      size_t strpid = i++;
      if( self->state[strpid] != block_t::Empty ) continue;
      writer->Read( self->blkid, strpid, self->stripes[strpid],
                         RepairTool::update_callback( self, writer, strpid ) ,0 , false);
      self->state[strpid] = block_t::Loading;
      ++loadingcnt;
    }

    return true;
}

//-----------------------------------------------------------------------
// Get a callback for read operation
//-----------------------------------------------------------------------
callback_t RepairTool::update_callback(std::shared_ptr<block_t> &self, RepairTool *tool,
		size_t strpid) {
	return [self, tool, strpid](const XrdCl::XRootDStatus &st, const uint32_t &length) mutable {
		std::unique_lock<std::mutex> lck(tool->blkmtx);
		self->state[strpid] = st.IsOK() ? self->Valid : self->Missing;
		if(st.IsOK()){
			tool->block->stripes[strpid].resize(length);
		}
		//------------------------------------------------------------
		// Check if we need to do any error correction (either for
		// the current stripe, or any other stripe)
		//------------------------------------------------------------
		if(!error_correction(self, tool)){
			tool->currentBlockChecked++;
			lck.unlock();
			tool->CheckBlock();
		}
	};
}

//-----------------------------------------------------------------------
// Get a callback for read operation
//-----------------------------------------------------------------------
callback_t RepairTool::read_callback(std::shared_ptr<ThreadEndSemaphore> sem, size_t blkid, size_t strpid, RepairTool *tool) {
	return [tool, sem, blkid, strpid](const XrdCl::XRootDStatus &st, const uint32_t &length) mutable {
		if(sem != nullptr && !st.IsOK()){
			std::cout << "Corruption in block " << blkid << " and stripe " << strpid << "\nFile: "
					<< tool->objcfg.GetFileName(blkid, strpid) << "\n" << std::flush;
			tool->st->status = XrdCl::stError;
		}
	};
}

void RepairTool::CheckFile(XrdCl::ResponseHandler *handler){
	XrdCl::SyncResponseHandler handler1;
	TryOpen(&handler1, XrdCl::OpenFlags::Read, 0);
	handler1.WaitForResponse();
	st = handler1.GetStatus();
	if (handler1.GetStatus()->IsOK())
	{

		auto itr = redirectionMap.begin();
		for (; itr != redirectionMap.end(); ++itr)
		{
			const std::string &oldUrl = itr->first;
			std::cout << "Archive contains some damaged metadata: " << oldUrl << "\n" << std::flush;
			InvalidateReplaceArchive(oldUrl, readDataarchs[oldUrl]);
			st->status = XrdCl::stError;
		}
	}
	// do the read for each strpid and blkid but with different callback func
	uint64_t numBlocks = ceil((filesize / (float) objcfg.chunksize) / objcfg.nbdata);
	std::vector<std::shared_ptr<buffer_t>> buffers;
	auto sem = std::make_shared<XrdSysSemaphore>(0);
	{
		// switching out of this context will remove the own reference to ptr
		std::shared_ptr<ThreadEndSemaphore> ptr = std::make_shared<
				ThreadEndSemaphore>(sem);

		for (size_t blkid = 0; blkid < numBlocks; blkid++)
		{
			for (size_t strpid = 0; strpid < objcfg.nbdata + objcfg.nbparity;
					strpid++)
			{
				std::shared_ptr<buffer_t> buffer = std::make_shared<buffer_t>();
				buffer->reserve(objcfg.chunksize);
				buffers.push_back(buffer);
				Read(blkid, strpid, *buffer,
						RepairTool::read_callback(ptr, blkid, strpid, this));
			}
		}
	}
	sem->Wait();
	for(size_t u = 0; u < objcfg.plgr.size(); u++){
		writeDataarchs[objcfg.GetDataUrl(u)]= readDataarchs[objcfg.GetDataUrl(u)];
	}
	XrdCl::SyncResponseHandler handler2;
	CloseAllArchives(&handler2);
	handler2.WaitForResponse();
	if (handler)
	{
		handler->HandleResponse(new XrdCl::XRootDStatus(*st), nullptr);
	}
}

void RepairTool::RepairFile(bool checkAgainAfterRepair, XrdCl::ResponseHandler *handler) {
	std::cout<<"Repair called with " << (int)objcfg.nbchunks << " chunks of which " << (int)objcfg.nbdata << "are data\n"<<std::flush;

	XrdCl::SyncResponseHandler handler1;
	RepairTool::OpenInUpdateMode(&handler1);
	handler1.WaitForResponse();

	std::cout<<"Opened archives\n"<<std::flush;
	/*for(auto it = urlmap.begin(); it != urlmap.end(); it++){
		std::cout << "Map key " << (*it).first << " to " << (*it).second << "\n" << std::flush;
	}*/
	checkAfterRepair = checkAgainAfterRepair;

	chunksRepaired.store(0);
	currentBlockChecked = 0;
	repairFailed = false;

	std::unique_lock<std::mutex> lk(finishedRepairMutex);
	CheckBlock();
	repairVar.wait(lk, [this]{return this->finishedRepair && this->chunkRepairsWritten.load() == this->chunksRepaired.load();});
	lk.unlock();

	XrdCl::SyncResponseHandler handlerClose;
	CloseAllArchives(&handlerClose, 0);
	handlerClose.WaitForResponse();

	std::cout<<"Archives closed\n"<< std::flush;


		for(size_t i = 0; i < objcfg.nbchunks; i++){
			if(redirectionMap.find(objcfg.GetDataUrl(i)) != redirectionMap.end()){
				objcfg.plgr[i] = std::string(redirectionMap[objcfg.GetDataUrl(i)]);
			}
		}
		for(auto it = objcfg.plgr.begin(); it != objcfg.plgr.end(); it++){
					std::cout << "Host at " << *it << "\n" << std::flush;
				}
		for(auto it = redirectionMap.begin(); it != redirectionMap.end(); it++)
							std::cout << "Redirect from " << it->first << " to " << it->second << "\n" << std::flush;


		std::cout << "\nCHUNKS REPAIRED: " << chunksRepaired.load() << "\n" << std::flush;
		if(repairFailed){
			std::cout << "Repair failed at some point!\n"<<std::flush;
		}

}

void RepairTool::CheckBlock() {
	std::unique_lock<std::mutex> lck(blkmtx);
	redirectMapOffset = 0;
	size_t totalBlocks = 0;
	if (objcfg.nomtfile) {
		totalBlocks = lstblk+1;
	} else {
		totalBlocks = lstblk + 1;
	}
	size_t blkid = currentBlockChecked;
	if (blkid < totalBlocks) {
		if (!block || block->blkid != blkid)
			block = std::make_shared<block_t>(blkid, reader, objcfg);
		auto blk = block;
		// initiates the read for each stripe, at some point CheckBlock will be called again from somewhere else.
		if (!error_correction(blk, this)) {
			currentBlockChecked++;
			lck.unlock();
			CheckBlock();
		}
	}
	else {
		std::unique_lock<std::mutex> lk(finishedRepairMutex);
		std::cout << "Finished Repair init\n"<<std::flush;
		finishedRepair = true;
		repairVar.notify_all();
	}

}

XrdCl::XRootDStatus RepairTool::WriteChunk(std::shared_ptr<block_t> blk, size_t strpid){
	auto blkid = blk->blkid;
	std::string fn = objcfg.GetFileName( blkid, strpid );

	std::string url;
	// if the block/stripe does not exist it means we are reading passed the end of the file
	auto itr = urlmap.find(fn);
	if (itr == urlmap.end())
	{
		std::cout << "Couldn't locate file " << fn << " for writing\n" << std::flush;
		auto it = redirectionMap.begin();
		int u = 0;
		while(it != redirectionMap.end()){
			if(u >= redirectMapOffset) break;
			u++;
			it++;
		}
		url = it->first;
		redirectMapOffset++;
		std::cout << "Replacing it with host " << url << "\n" << std::flush;
	}
	else url = itr->second;
	if(writeDataarchs[url]->IsOpen()){
		std::string fn = objcfg.GetFileName( blkid, strpid );
		auto it = writeDataarchs[url]->cdmap.find(fn);
		uint64_t offset = 0;
		int32_t actualSize = objcfg.chunksize;
		if (strpid < objcfg.nbdata && filesize < blkid * objcfg.chunksize * objcfg.nbdata
						+ (strpid + 1) * objcfg.chunksize)
		{
			actualSize = filesize - (blkid * objcfg.chunksize * objcfg.nbdata
							+ strpid * objcfg.chunksize);
			if (actualSize < 0)
				actualSize = 0;
		}
		// if we are in a parity stripe and only the first data stripe has size > 0, the parity stripe has that same size ( < chunksize)
		if(strpid >= objcfg.nbdata && filesize < blkid * objcfg.chunksize * objcfg.nbdata + objcfg.chunksize){
			actualSize = filesize - (blkid * objcfg.chunksize * objcfg.nbdata);
			if (actualSize < 0)
				actualSize = 0;
		}
		std::cout << "Actual Write Size: " << actualSize
				<< " because file size is " << filesize << "\n" << std::flush;

		auto pipehndl = [=](const XrdCl::XRootDStatus &st) {
			// increase the written counter by one (atomic)
			if(!st.IsOK()){
				std::cout << "Write to " << url << " failed: " << st.code << "\n" << std::flush;
			}
			std::unique_lock<std::mutex> lk(finishedRepairMutex);
			chunkRepairsWritten.fetch_add(1);

			std::cout << "Write pipeline ended "<<chunkRepairsWritten.load() << "\n" << std::flush;
			repairVar.notify_all();
			    	};

		if(it != writeDataarchs[url]->cdmap.end()){

			// the file exists, so we overwrite it (but make a copy of the current state for testing purposes)
			uint32_t chksum = objcfg.digest(0, &(blk->stripes[strpid][0]), actualSize);
			std::stringstream ss;
			ss << "cp " << url << " " << url<< strpid;
			std::string s = ss.str();
			system(s.data());
			//return writeDataarchs[url]->WriteFileInto(fn, offset, actualSize, chksum, &(blk->stripes[strpid][0]), nullptr, 0);
			XrdCl::Pipeline p = XrdCl::WriteIntoFile(XrdCl::Ctx<XrdCl::ZipArchive>(*writeDataarchs[url]),
					fn, offset, actualSize, chksum, &(blk->stripes[strpid][0]), 0)
							| XrdCl::Final(pipehndl);
			XrdCl::Async(std::move(p));
			return XrdCl::XRootDStatus();
		}
		else{

			uint32_t chksum = objcfg.digest(0, &(blk->stripes[strpid][0]), actualSize);
			// the file doesnt exist at all, so we append it to the archive (likely a completely new archive)
			XrdCl::Pipeline p = XrdCl::AppendFile(
					XrdCl::Ctx<XrdCl::ZipArchive>(*writeDataarchs[url]),
					fn, chksum, actualSize, &(blk->stripes[strpid][0]))
				| XrdCl::Final(pipehndl);
			XrdCl::Async(std::move(p));
			return XrdCl::XRootDStatus();
		}
	}
	std::unique_lock<std::mutex> lk(finishedRepairMutex);
				std::cout << "Write call failed\n" << std::flush;
				chunkRepairsWritten.fetch_add(1);
				repairVar.notify_all();
	return XrdCl::XRootDStatus(XrdCl::stError, "Couldnt write to archive");
}

//---------------------------------------------------------------------------
// Get a buffer with metadata (CDFH and EOCD records)
//---------------------------------------------------------------------------
XrdZip::buffer_t RepairTool::GetMetadataBuffer()
{
  using namespace XrdZip;

  const size_t cdcnt = objcfg.plgr.size();
  std::vector<buffer_t> buffs; buffs.reserve( cdcnt ); // buffers with raw data
  std::vector<LFH> lfhs; lfhs.reserve( cdcnt );        // LFH records
  std::vector<CDFH> cdfhs; cdfhs.reserve( cdcnt );     // CDFH records

  //-------------------------------------------------------------------------
  // prepare data structures (LFH and CDFH records)
  //-------------------------------------------------------------------------
  uint64_t offset = 0;
  uint64_t cdsize = 0;
  mode_t mode = S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH;
  for( size_t i = 0; i < cdcnt; ++i )
  {
    std::string fn = std::to_string( i );                          // file name (URL of the data archive)
    buffer_t buff( writeDataarchs[objcfg.GetDataUrl(i)]->GetCD() );                        // raw data buffer (central directory of the data archive)
    uint32_t cksum = objcfg.digest( 0, buff.data(), buff.size() ); // digest (crc) of the buffer
    lfhs.emplace_back( fn, cksum, buff.size(), time( 0 ) );        // LFH record for the buffer
    LFH &lfh = lfhs.back();
    cdfhs.emplace_back( &lfh, mode, offset );                      // CDFH record for the buffer
    offset += LFH::lfhBaseSize + fn.size() + buff.size();          // shift the offset
    cdsize += cdfhs.back().cdfhSize;                               // update central directory size
    buffs.emplace_back( std::move( buff ) );                       // keep the buffer for later
  }

  uint64_t zipsize = offset + cdsize + EOCD::eocdBaseSize;
  buffer_t zipbuff; zipbuff.reserve( zipsize );

  //-------------------------------------------------------------------------
  // write into the final buffer LFH records + raw data
  //-------------------------------------------------------------------------
  for( size_t i = 0; i < cdcnt; ++i )
  {
    lfhs[i].Serialize( zipbuff );
    std::copy( buffs[i].begin(), buffs[i].end(), std::back_inserter( zipbuff ) );
  }
  //-------------------------------------------------------------------------
  // write into the final buffer CDFH records
  //-------------------------------------------------------------------------
  for( size_t i = 0; i < cdcnt; ++i )
    cdfhs[i].Serialize( zipbuff );
  //-------------------------------------------------------------------------
  // prepare and write into the final buffer the EOCD record
  //-------------------------------------------------------------------------
  EOCD eocd( offset, cdcnt, cdsize );
  eocd.Serialize( zipbuff );

  return zipbuff;
}

void RepairTool::CloseAllArchives(XrdCl::ResponseHandler *handler, uint16_t timeout){
	//-------------------------------------------------------------------------
	    // First, check the global status, if we are in an error state just
	    // fail the request.
	    //-------------------------------------------------------------------------

	    const size_t size = objcfg.plgr.size();
	    //-------------------------------------------------------------------------
	    // prepare the metadata (the Central Directory of each data ZIP)
	    //-------------------------------------------------------------------------
	    auto zipbuff = objcfg.nomtfile ? nullptr :
	                   std::make_shared<XrdZip::buffer_t>( GetMetadataBuffer() );
	    //-------------------------------------------------------------------------
	    // prepare the pipelines ...
	    //-------------------------------------------------------------------------
	    std::vector<XrdCl::Pipeline> closes;
	    std::vector<XrdCl::Pipeline> save_metadata;
	    closes.reserve( size );
	    std::string closeTime = std::to_string( time(NULL) );

	    // since the filesize was already determined in previous writes we keep it
	    //XrdCl::Pipeline p2 = ReadSize(0);
	    //XrdCl::Async(std::move(p2));

	    std::vector<XrdCl::xattr_t> xav{ {"xrdec.filesize", std::to_string(filesize)},
	                                     {"xrdec.strpver", closeTime.c_str()},
										 {"xrdec.corrupted", std::to_string(0)}};

	    for( size_t i = 0; i < size; ++i )
	    {
	    // we only close written archives since any replaced archives are closed in InvalidateReplaceArchive method
	      //-----------------------------------------------------------------------
	      // close ZIP archives with data
	      //-----------------------------------------------------------------------
	      if( writeDataarchs[objcfg.GetDataUrl(i)]->IsOpen() )
	      {
	    	  XrdCl::Pipeline p = XrdCl::SetXAttr( writeDataarchs[objcfg.GetDataUrl(i)]->GetFile(), xav )
	                          | XrdCl::CloseArchive( *writeDataarchs[objcfg.GetDataUrl(i)] );
	        closes.emplace_back( std::move( p ) );
	      }
	      else{
	    	  std::cout<< "Archive " << objcfg.GetDataUrl(i) << " not open but rather " << writeDataarchs[objcfg.GetDataUrl(i)]->openstage << "\n" << std::flush;
	    	  XrdCl::Pipeline p = XrdCl::SetXAttr( writeDataarchs[objcfg.GetDataUrl(i)]->GetFile(), xav )
	    	  	                          | XrdCl::CloseArchive( *writeDataarchs[objcfg.GetDataUrl(i)] );
	    	  closes.emplace_back( std::move( p ) );
	      }
	      //-----------------------------------------------------------------------
	      // replicate the metadata
	      //-----------------------------------------------------------------------
	      if( zipbuff )
	      {
	        std::string url = objcfg.GetMetadataUrl( i );
	        metadataarchs.emplace_back( std::make_shared<XrdCl::File>(
	            Config::Instance().enable_plugins ) );
	        XrdCl::Pipeline p = XrdCl::Open( *metadataarchs[i], url, XrdCl::OpenFlags::New | XrdCl::OpenFlags::Write )
	                          | XrdCl::Write( *metadataarchs[i], 0, zipbuff->size(), zipbuff->data() )
	                          | XrdCl::Close( *metadataarchs[i] )
	                          | XrdCl::Final( [zipbuff]( const XrdCl::XRootDStatus& ){ } );

	        save_metadata.emplace_back( std::move( p ) );
	      }
	    }

	    auto pipehndl = [=](const XrdCl::XRootDStatus &st) { // set the central directories in ZIP archives (if we use metadata files)
	    		auto itr = writeDataarchs.begin();
	    		for (; itr != writeDataarchs.end(); ++itr) {
	    			auto &zipptr = itr->second;
	    			auto url = itr->first;
	    			if(zipptr->openstage != XrdCl::ZipArchive::None){
	    				std::cout << "Archive wasn't properly closed: " << url << ", status " << zipptr->openstage << "\n" << std::flush;
	    				XrdCl::StatInfo* info = new XrdCl::StatInfo();
	    				XrdCl::XRootDStatus s = zipptr->archive.Stat(false, info);
	    				if(s.IsOK()){
	    					std::cout << "File state: " << info->GetFlags();
	    				}
	    				else std::cout << "Couldn't stat file\n"<<std::flush;
	    			}
	    		}
	    		if (handler)
	    			handler->HandleResponse(new XrdCl::XRootDStatus(st), nullptr);
	    	};

	    //-------------------------------------------------------------------------
	    // If we were instructed not to create the the additional metadata file
	    // do the simplified close
	    //-------------------------------------------------------------------------
	    if( save_metadata.empty() )
	    {
	      XrdCl::Pipeline p = XrdCl::Parallel( closes ).AtLeast( objcfg.nbchunks )| XrdCl::Final(pipehndl);
	      XrdCl::Async( std::move( p ), timeout );
	      return;
	    }

	    //-------------------------------------------------------------------------
	    // compose closes & save_metadata:
	    //  - closes must be successful at least for #data + #parity
	    //  - save_metadata must be successful at least for #parity + 1
	    //-------------------------------------------------------------------------
	    XrdCl::Pipeline p = XrdCl::Parallel(
	        XrdCl::Parallel( closes ).AtLeast( objcfg.nbchunks ),
	        XrdCl::Parallel( save_metadata ).AtLeast( objcfg.nbparity + 1 )
	      );
	    XrdCl::Async( std::move( p ), timeout );
}

void RepairTool::TryOpen(XrdCl::ResponseHandler *handler, XrdCl::OpenFlags::Flags flags, uint16_t timeout){
	const size_t size = objcfg.plgr.size();
	std::vector<XrdCl::Pipeline> opens;
	opens.reserve(size);
	std::vector<XrdCl::Pipeline> healthRead;
	healthRead.reserve(size);
	for (size_t i = 0; i < size; ++i)
	{
		// generate the URL
		std::string url = objcfg.GetDataUrl(i);
		auto archive = std::make_shared<XrdCl::ZipArchive>(
				Config::Instance().enable_plugins);
		// create the file object
		readDataarchs.emplace(url, archive);
		//writeDataarchs.emplace(url, archive);
		// open the archive
		if (objcfg.nomtfile)
		{
			opens.emplace_back(
					XrdCl::OpenArchive(*readDataarchs[url], url,
							flags));
		}
		else
			opens.emplace_back(OpenOnly(*readDataarchs[url], url, true));
		healthRead.emplace_back(CheckHealthExists(i));
	}

	auto pipehndl =
			[=](const XrdCl::XRootDStatus &st)
			{ // set the central directories in ZIP archives (if we use metadata files)
						auto itr = readDataarchs.begin();
						for (; itr != readDataarchs.end(); ++itr)
						{
							const std::string &url = itr->first;
							auto &zipptr = itr->second;
							if (zipptr->openstage == XrdCl::ZipArchive::NotParsed)
							zipptr->SetCD(metadata[url]);
							else if (zipptr->openstage != XrdCl::ZipArchive::Done)
							{
								ReplaceURL(url);
								if(!metadata.empty())
								AddMissing(metadata[url]);
							}
							auto itr = zipptr->cdmap.begin();
							for (; itr != zipptr->cdmap.end(); ++itr)
							{
								try
								{
									size_t blknb = fntoblk(itr->first);
									urlmap.emplace(itr->first, url);
									if (blknb > lstblk)
									lstblk = blknb;
								}
								catch (std::invalid_argument&)
								{
									std::cout << "Invalid file name detected\n" << std::flush;
								}
							}
						}
						metadata.clear();
						auto sem = std::make_shared<XrdSysSemaphore>(0);
						{
							std::shared_ptr<ThreadEndSemaphore> ptr = std::make_shared<ThreadEndSemaphore>(sem);
							// Check that all LFH and CDFH are correct
							CheckAllMetadata(ptr);
						}
						sem->Wait();
						// call user handler
						if (handler)
						handler->HandleResponse(new XrdCl::XRootDStatus(st), nullptr);
					};
	// in parallel open the data files and read the metadata
	XrdCl::Pipeline p =
			objcfg.nomtfile ?
					XrdCl::Parallel(opens).AtLeast(objcfg.nbdata)
							| XrdCl::Parallel(healthRead).AtLeast(0)
							| ReadSize(0) | XrdCl::Final(pipehndl) :
					XrdCl::Parallel(ReadMetadata(0),
							XrdCl::Parallel(opens).AtLeast(objcfg.nbdata))
							>> pipehndl;
	XrdCl::Async(std::move(p), timeout);
}

void RepairTool::OpenInUpdateMode(XrdCl::ResponseHandler *handler,
		uint16_t timeout) {
	XrdCl::SyncResponseHandler handler1;
	TryOpen(&handler1, XrdCl::OpenFlags::Update, 0);
	handler1.WaitForResponse();
	if (handler1.GetStatus()->IsOK())
	{
		auto readItr = readDataarchs.begin();
		for(; readItr != readDataarchs.end(); ++readItr){
			writeDataarchs.emplace(readItr->first, readItr->second);
		}
		auto itr = redirectionMap.begin();
		size_t index = 0;
		for (; itr != redirectionMap.end(); ++itr)
		{
			const std::string &oldUrl = itr->first;
			// the redirection map saves the plgr without the "test.txt" file name for easier future replacement
			const std::string &newUrl = objcfg.GetReplacementUrl(index);
			auto newArch = std::make_shared<XrdCl::ZipArchive>(
					Config::Instance().enable_plugins);
			newArch->OpenArchive(newUrl,
					XrdCl::OpenFlags::New | XrdCl::OpenFlags::Write, nullptr,
					0);
			writeDataarchs[oldUrl] = newArch;
			std::cout << "Creating new archive pointing from url " << oldUrl
					<< " to " << newUrl << "\n" << std::flush;

			InvalidateReplaceArchive(oldUrl, readDataarchs[oldUrl]);

			index++;
		}
	}
	if (handler)
	{
		handler->HandleResponse(handler1.GetStatus(), nullptr);
	}
}

XrdCl::Pipeline RepairTool::CheckHealthExists(size_t index){
	  std::string url = objcfg.GetDataUrl( index );
	  		  return XrdCl::ListXAttr(readDataarchs[url]->GetFile()) >>
	  				  [index, url, this] (XrdCl::XRootDStatus &st, std::vector<XrdCl::XAttr> attrs){
	  			  	  	 for(auto it = attrs.begin(); it != attrs.end(); it++){
	  			  	  		 if(it->name == "xrdec.corrupted"){
	  			  	  			 XrdCl::Pipeline::Replace(ReadHealth(index));
	  			  	  		 }
	  			  	  	 }
	  		  };
  }

  XrdCl::Pipeline RepairTool::ReadHealth(size_t index){
		  std::string url = objcfg.GetDataUrl( index );
		  return XrdCl::GetXAttr( readDataarchs[url]->GetFile(), "xrdec.corrupted" ) >>
	          [index, url, this]( XrdCl::XRootDStatus &st, std::string &damage)
	          {
			  //std::cout << "Attempt to read health\n" << std::flush;
				if (st.IsOK()) {
					try {
						int damaged = std::stoi(damage);
						if (damaged > 0)
							this->readDataarchs[url]->openstage = XrdCl::ZipArchive::Error;
					} catch (std::invalid_argument&) {
						return;
					}
				}
				return;
	          };
  }

void RepairTool::CheckAllMetadata(std::shared_ptr<ThreadEndSemaphore> sem) {
	uint64_t numBlocks = ceil(
			(filesize / (float)objcfg.chunksize) / objcfg.nbdata);
	for (size_t blkid = 0; blkid < numBlocks; blkid++) {
		std::vector<size_t> stripesMissing;
		for (size_t strpid = 0; strpid < objcfg.nbdata + objcfg.nbparity;
				strpid++) {
			std::string fn = objcfg.GetFileName(blkid, strpid);
				// if the block/stripe does not exist it means we are reading passed the end of the file
				auto itr = urlmap.find(fn);
				if (itr == urlmap.end()) {
					stripesMissing.push_back(strpid);
				}
				else CompareLFHToCDFH(sem, blkid, strpid);
		}
		if(stripesMissing.size() > 0){
			// replace all URLs that did not receive a stripe of this block (since they must be damaged)
			// if they are already in the redirection map, don't do anything.
			std::vector<std::string> dataarchUrls;
			for(auto it = readDataarchs.begin(); it != readDataarchs.end(); it++){
				bool missingEntry = true;
				for (size_t strpid = 0; strpid < objcfg.nbdata + objcfg.nbparity;
									strpid++) {
								std::string fn = objcfg.GetFileName(blkid, strpid);
								auto itr = urlmap.find(fn);
								if (itr != urlmap.end()&&itr->second == it->first) {
									missingEntry = false;
									break;
								}
							}
				if(missingEntry)
					dataarchUrls.push_back(it->first);
			}
			for(auto it = dataarchUrls.begin(); it != dataarchUrls.end(); it++){
				if(redirectionMap.find(*it) == redirectionMap.end()){
					ReplaceURL(*it);
				}
			}
		}
	}
}

void RepairTool::InvalidateReplaceArchive(const std::string &url, std::shared_ptr<XrdCl::ZipArchive> zipptr){
	XrdCl::Pipeline p;
	if(zipptr->IsOpen()){
		std::cout << "Closing and invalidating " << url << "\n" << std::flush;
		std::vector<XrdCl::xattr_t> xav{ {"xrdec.corrupted", std::to_string(1)} };
		p = XrdCl::SetXAttr( zipptr->archive, xav )
		                          | XrdCl::CloseArchive( *zipptr);
		XrdCl::Async(std::move(p), 0);
	}

}

void RepairTool::CompareLFHToCDFH(std::shared_ptr<ThreadEndSemaphore> sem, uint16_t blkid, uint16_t strpid){
	std::string fn = objcfg.GetFileName(blkid, strpid);
	// if the block/stripe does not exist it means we are reading passed the end of the file
	auto itr = urlmap.find(fn);
	if (itr == urlmap.end()) {
		return;
	}
	// get the URL of the ZIP archive with the respective data
	const std::string &url = itr->second;
	if (redirectionMap.find(url) != redirectionMap.end()) {
		// the url was already marked as non existant / not reachable, skip
		return;
	}
	if(readDataarchs.find(url) == readDataarchs.end()){
		// the url doesn't exist / not reachable, we need to redirect.
		// Should be caught by previous OpenInUpdateMode
		ReplaceURL(url);
	}
	std::shared_ptr<XrdCl::ZipArchive> &zipptr = readDataarchs[url];
	if(zipptr->openstage != XrdCl::ZipArchive::Done){
		ReplaceURL(url);
	}
	auto cditr = zipptr->cdmap.find(fn);
	if (cditr == zipptr->cdmap.end()) {
		ReplaceURL(url);
		return;
	}

	XrdZip::CDFH *cdfh = zipptr->cdvec[cditr->second].get();
	uint32_t offset = XrdZip::CDFH::GetOffset(*cdfh);
	uint64_t cdOffset =
			zipptr->zip64eocd ?
					zipptr->zip64eocd->cdOffset : zipptr->eocd->cdOffset;
	uint64_t nextRecordOffset =
			(cditr->second + 1 < zipptr->cdvec.size()) ?
					XrdZip::CDFH::GetOffset(*zipptr->cdvec[cditr->second + 1]) :
					cdOffset;
	auto readSize = (nextRecordOffset - offset) - cdfh->uncompressedSize;
	//- objcfg.chunksize;
	std::shared_ptr<buffer_t> lfhbuf;
	lfhbuf = std::make_shared<buffer_t>();
	lfhbuf->reserve(readSize);

	std::shared_ptr<ThreadEndSemaphore> localSem(sem);
	XrdCl::Pipeline p = XrdCl::Read(zipptr->archive, offset, readSize,
			lfhbuf->data(), 0)
			>> [=](XrdCl::XRootDStatus &st, XrdCl::ChunkInfo &ch) {
				if (st.IsOK() && localSem != nullptr) {
					try {
						XrdZip::LFH lfh(lfhbuf->data(), readSize);
						if (lfh.ZCRC32 != cdfh->ZCRC32|| lfh.compressedSize != cdfh->compressedSize
								|| lfh.compressionMethod != cdfh->compressionMethod
								|| lfh.extraLength != cdfh->extraLength
								|| lfh.filename != cdfh->filename
								|| lfh.filenameLength != cdfh->filenameLength
								|| lfh.generalBitFlag != cdfh->generalBitFlag
								|| lfh.minZipVersion != cdfh->minZipVersion
								|| lfh.uncompressedSize != cdfh->uncompressedSize) {
							// metadata damaged, mark archive as damaged and replace url
							ReplaceURL(url);
							return;
						}
						//---------------------------------------------------
						// All is good
						//---------------------------------------------------
						else {
							//std::cout << "Header checked " << fn << "\n" << std::flush;
							return;
						}
					} catch (const std::exception &e) {
						// Couldn't parse metadata, metadata damaged, same case.
						ReplaceURL(url);
						return;
					}
				} else {
					// Couldn't even read from file?! Shouldn't happen
					ReplaceURL(url);
					return;
				}
			};
	Async( std::move( p ), 0 );

}

void RepairTool::ReplaceURL(const std::string &url){
	urlMutex.lock();
	if(redirectionMap.find(url) != redirectionMap.end()){
		urlMutex.unlock();
		return;
	}
	if (objcfg.plgrReplace.size() > currentReplaceIndex) {
			const std::string replacementPlgr = objcfg.plgrReplace[currentReplaceIndex];
			// save a mapping from old to new urls so user knows where actual data is
			redirectionMap[url] = replacementPlgr;
			currentReplaceIndex++;
		}else{
			// TODO: throw error if there's no new url
			std::cout << "Critical error, can't find replacement host for " << url << "\n" << std::flush;
			redirectionMap[url] = "null";
		}
	urlMutex.unlock();
}

//-----------------------------------------------------------------------
// Add all the entries from given Central Directory to missing
//-----------------------------------------------------------------------
void RepairTool::AddMissing(const buffer_t &cdbuff) {
	const char *buff = cdbuff.data();
	size_t size = cdbuff.size();
	// parse Central Directory records
	XrdZip::cdvec_t cdvec;
	XrdZip::cdmap_t cdmap;
	std::tie(cdvec, cdmap) = XrdZip::CDFH::Parse(buff, size);
	auto itr = cdvec.begin();
	for (; itr != cdvec.end(); ++itr) {
		XrdZip::CDFH &cdfh = **itr;
		missing.insert(cdfh.filename);
	}
}

//-----------------------------------------------------------------------
//! Read size from xattr
//!
//! @param index : placement's index
//-----------------------------------------------------------------------
XrdCl::Pipeline RepairTool::ReadSize(size_t index) {
	std::string url = objcfg.GetDataUrl(index);
	return XrdCl::GetXAttr(readDataarchs[url]->GetFile(), "xrdec.filesize")
			>> [index, this](XrdCl::XRootDStatus &st, std::string &size) {
				if (!st.IsOK()) {
					//-------------------------------------------------------------
					// Check if we can recover the error or a diffrent location
					//-------------------------------------------------------------
					if (index + 1 < objcfg.plgr.size())
						XrdCl::Pipeline::Replace(ReadSize(index + 1));
					return;
				}
				auto newFileSize = std::stoull(size);
				if(newFileSize == 0){
					if (index + 1 < objcfg.plgr.size())
											XrdCl::Pipeline::Replace(ReadSize(index + 1));
										return;
				}
				filesize = newFileSize;
			};
}

//-----------------------------------------------------------------------
// Parse metadata from chunk info object
//-----------------------------------------------------------------------
bool RepairTool::ParseMetadata(XrdCl::ChunkInfo &ch) {
	const size_t mincnt = objcfg.nbdata + objcfg.nbparity;
	const size_t maxcnt = objcfg.plgr.size();

	char *buffer = reinterpret_cast<char*>(ch.buffer);
	size_t length = ch.length;

	for (size_t i = 0; i < maxcnt; ++i) {
		uint32_t signature = XrdZip::to<uint32_t>(buffer);
		if (signature != XrdZip::LFH::lfhSign) {
			if (i + 1 < mincnt)
				return false;
			break;
		}
		XrdZip::LFH lfh(buffer);
		// check if we are not reading passed the end of the buffer
		if (lfh.lfhSize + lfh.uncompressedSize > length)
			return false;
		buffer += lfh.lfhSize;
		length -= lfh.lfhSize;
		// verify the checksum
		uint32_t crc32val = objcfg.digest(0, buffer, lfh.uncompressedSize);
		if (crc32val != lfh.ZCRC32)
			return false;
		// keep the metadata
		std::string url = objcfg.GetDataUrl(std::stoull(lfh.filename));
		metadata.emplace(url, buffer_t(buffer, buffer + lfh.uncompressedSize));
		buffer += lfh.uncompressedSize;
		length -= lfh.uncompressedSize;
	}

	return true;
}

//-----------------------------------------------------------------------
// Read metadata for the object
//-----------------------------------------------------------------------
XrdCl::Pipeline RepairTool::ReadMetadata(size_t index) {
	const size_t size = objcfg.plgr.size();
	// create the File object
	auto file = std::make_shared<XrdCl::File>(
			Config::Instance().enable_plugins);
	// prepare the URL for Open operation
	std::string url = objcfg.GetMetadataUrl(index);
	// arguments for the Read operation
	XrdCl::Fwd<uint32_t> rdsize;
	XrdCl::Fwd<void*> rdbuff;

	return XrdCl::Open(*file, url, XrdCl::OpenFlags::Read)
			>> [=](XrdCl::XRootDStatus &st, XrdCl::StatInfo &info) mutable {
				if (!st.IsOK()) {
					if (index + 1 < size)
						XrdCl::Pipeline::Replace(ReadMetadata(index + 1));
					return;
				}
				// prepare the args for the subsequent operation
				rdsize = info.GetSize();
				rdbuff = new char[info.GetSize()];
			}
			| XrdCl::Read(*file, 0, rdsize, rdbuff)
					>> [=](XrdCl::XRootDStatus &st, XrdCl::ChunkInfo &ch) {
						if (!st.IsOK()) {
							if (index + 1 < size)
								XrdCl::Pipeline::Replace(
										ReadMetadata(index + 1));
							return;
						}
						// now parse the metadata
						if (!ParseMetadata(ch)) {
							if (index + 1 < size)
								XrdCl::Pipeline::Replace(
										ReadMetadata(index + 1));
							return;
						}
					} | XrdCl::Close(*file) >> [](XrdCl::XRootDStatus &st) {
				if (!st.IsOK())
					XrdCl::Pipeline::Ignore(); // ignore errors, we don't really care
			} | XrdCl::Final([rdbuff, file](const XrdCl::XRootDStatus&) {
				// deallocate the buffer if necessary
				if (rdbuff.Valid()) {
					char *buffer = reinterpret_cast<char*>(*rdbuff);
					delete[] buffer;
				}
			});
}

void RepairTool::Read( size_t blknb, size_t strpnb, buffer_t &buffer, callback_t cb, uint16_t timeout, bool exactControl )
{
  // generate the file name (blknb/strpnb)
  std::string fn = objcfg.GetFileName( blknb, strpnb );
  // if the block/stripe does not exist it means we are reading passed the end of the file
  auto itr = urlmap.find( fn );
  if( itr == urlmap.end() )
  {
    auto st = !IsMissing( fn ) ? XrdCl::XRootDStatus() :
              XrdCl::XRootDStatus( XrdCl::stError, XrdCl::errNotFound );
    ThreadPool::Instance().Execute( cb, st, 0 );
    return;
  }
  // get the URL of the ZIP archive with the respective data
  const std::string &url = itr->second;

  if(redirectionMap.find(url) != redirectionMap.end()){
	  auto st = XrdCl::XRootDStatus(XrdCl::stError, XrdCl::errRedirect);
	  ThreadPool::Instance().Execute(cb, st, 0);
	  return;
  }

  // get the ZipArchive object
  auto &zipptr = readDataarchs[url];
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
  XrdCl::Async( XrdCl::ReadFrom( *zipptr, fn, 0, rdsize, buffer.data() ) >>
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
																//- objcfg.chunksize;
							std::shared_ptr<buffer_t> lfhbuf;
							lfhbuf = std::make_shared<buffer_t>();
							lfhbuf->reserve(readSize);
							XrdCl::Pipeline p = XrdCl::Read(zipptr->archive, offset, readSize, lfhbuf->data(), timeout) >>
							                   [=]( XrdCl::XRootDStatus &st, XrdCl::ChunkInfo &ch )
							                   {
												if (st.IsOK()) {
													try {
														XrdZip::LFH lfh(
																lfhbuf->data());
														if (lfh.ZCRC32
																!= orgcksum){
															cb(
																	XrdCl::XRootDStatus(
																			XrdCl::stError,
																			"CRC LFH != CRC of CDFH"),
																	0);
															return;
														}
																if(lfh.compressedSize
																		!= cdfh->compressedSize){
																	cb(
																			XrdCl::XRootDStatus(
																					XrdCl::stError,
																					"CompSize LFH !=CDFH"),
																			0);
																	return;
																}
																if( lfh.compressionMethod
																		!= cdfh->compressionMethod
																|| lfh.extraLength
																		!= cdfh->extraLength
																|| lfh.filename
																		!= cdfh->filename
																|| lfh.filenameLength
																		!= cdfh->filenameLength
																|| lfh.generalBitFlag
																		!= cdfh->generalBitFlag
																|| lfh.minZipVersion
																		!= cdfh->minZipVersion){
																	cb(
																			XrdCl::XRootDStatus(
																					XrdCl::stError,
																					"LFH!=CDFH"),
																			0);
																	return;
																}
																if( lfh.uncompressedSize
																		!= cdfh->uncompressedSize) {
															cb(
																	XrdCl::XRootDStatus(
																			XrdCl::stError,
																			"uncompSize LFH !=CDFH"),
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
                  }, timeout );
}

//-----------------------------------------------------------------------
//! Check if chunk file name is missing
//-----------------------------------------------------------------------
bool RepairTool::IsMissing( const std::string &fn )
{
  // if the chunk is in the missing set return true
  if( missing.count( fn ) ) return true;
  // if we don't have a metadata file and the chunk exceeds last chunk
  // also return true
  try{
      if( objcfg.nomtfile && fntoblk( fn ) <= lstblk ) return true;}
      catch(...){}// otherwise return false
  return false;
}


}
