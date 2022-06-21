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
	//std::cout<<"Error Correction called with " << missingcnt << " missing and " << validcnt << " valid chunks\n"<<std::flush;
    if(validcnt == writer->objcfg.nbchunks){
    	std::cout<<"All Errors corrected\n"<<std::flush;
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
      return false;
    }
    //---------------------------------------------------------------------
    // Check if we can do the recovery right away
    //---------------------------------------------------------------------
    if( validcnt >= self->objcfg.nbdata && missingcnt > 0 )
    {
    	std::cout<<"Recovery possible, attempt starting\n"<<std::flush;

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
        std::cout<<"Write Chunk to disk: " << strpid << "\n"<<std::flush;
        std::cout << "Stripe " << (int)strpid << " with data ";
            		for(int i = 0; i < (int)writer->objcfg.chunksize; i++){
            			std::cout << self->stripes[strpid][i];
            		}
        writer->WriteChunk(self, strpid);
        if(writer->checkAfterRepair){
        	writer->Read( self->blkid, strpid, self->stripes[strpid],
        	                   RepairTool::update_callback( self, writer, strpid ) );
        }
        else{
        	writer->chunksRepaired++;
        	self->state[strpid] = block_t::Valid;
        	validcnt++;
        	if(validcnt == writer->objcfg.nbchunks){
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
      //std::cout<<"Load chunk of stripe " << strpid << " from disk\n"<<std::flush;
      writer->Read( self->blkid, strpid, self->stripes[strpid],
                         RepairTool::update_callback( self, writer, strpid ) );
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
	return [self, tool, strpid](const XrdCl::XRootDStatus &st, uint32_t) mutable {
		std::unique_lock<std::mutex> lck(tool->blkmtx);
		//std::cout << "Update callback for stripe " << strpid << " was " << (st.IsOK()?"successful":"unsuccessful")<<"\n" << std::flush;
		self->state[strpid] = st.IsOK() ? self->Valid : self->Missing;
		//------------------------------------------------------------
		// Check if we need to do any error correction (either for
		// the current stripe, or any other stripe)
		//------------------------------------------------------------
		if(!error_correction(self, tool)){
			tool->currentBlockChecked++;
			lck.unlock();
			tool->CheckBlock();
			std::cout<<"Out of Check Block of update_callback method\n"<< std::flush;
		}
	};
}

// might want to add a buffer so we don't create one every time?
/*void update_block_t(std::shared_ptr<block_t> &self, size_t strpid,
		uint32_t size,
		//callback_t                usrcb,
		uint16_t timeout) {
	std::unique_lock<std::mutex> lck(self->mtx);

	//---------------------------------------------------------------------
	// The cache is empty, we need to load the data
	//---------------------------------------------------------------------
	if (self->state[strpid] == self->Empty) {
		self->reader.Read(self->blkid, strpid, self->stripes[strpid],
				RepairTool::update_callback(self, strpid), timeout);
		self->state[strpid] = self->Loading;
	}
	//---------------------------------------------------------------------
	// The stripe is either corrupted or unreachable
	//---------------------------------------------------------------------
	if (self->state[strpid] == self->Missing) {
		if (!error_correction(self)) {
			//-----------------------------------------------------------------
			// Recovery was not possible, notify the user of the error
			//-----------------------------------------------------------------
			//usrcb(XrdCl::XRootDStatus(XrdCl::stError, XrdCl::errDataError), 0);
			return;
		}
		//-------------------------------------------------------------------
		// we fall through to the following if-statements that will handle
		// Recovering / Valid state
		//-------------------------------------------------------------------
	}
	//---------------------------------------------------------------------
	// The cache is loading or recovering, we don't have the data yet
	//---------------------------------------------------------------------
	if (self->state[strpid] == self->Loading || self->state[strpid] == self->Recovering) {
		return;
	}
	//---------------------------------------------------------------------
	// We do have the data so the chunk does not get written to disk again
	//---------------------------------------------------------------------
	if (self->state[strpid] == self->Valid) {
		//if (offset + size > self->stripes[strpid].size())
		//	size = self->stripes[strpid].size() - offset;
		//memcpy(usrbuff, self->stripes[strpid].data() + offset, size);
		//usrcb(XrdCl::XRootDStatus(), size);
		return false;
	}
	//---------------------------------------------------------------------
	// In principle we should never end up here, nevertheless if this
	// happens it is clearly an error ...
	//---------------------------------------------------------------------
	usrcb(XrdCl::XRootDStatus(XrdCl::stError, XrdCl::errInvalidOp), 0);
}*/

std::unique_ptr<ObjCfg> RepairTool::RepairFile(bool checkAgainAfterRepair, XrdCl::ResponseHandler *handler) {
	// pseudocode:
	// open file in update mode
	// create redirection table of hosts (index is stripe index)
	// for stripe i in stripes
	// check if host is reachable
	// if yes, insert host into table
	// if no, insert spare host into table

	// for block b in blocks (parallelized)
	// 	for chunk c in chunks
	//	check c
	//  if damaged/missing:
	// 		reconstruct c in b
	// 		cache b
	//		write c to host
	//		if re-read option:
	//		read c again and jump to "check c"
	// close file and return
	std::cout<<"Repair called with " << (int)objcfg.nbchunks << " chunks of which " << (int)objcfg.nbdata << "are data\n"<<std::flush;
	XrdCl::SyncResponseHandler handler1;
	RepairTool::OpenInUpdateMode(&handler1);
	handler1.WaitForResponse();



	std::cout<<"Open called\n"<<std::flush;
	for(auto it = urlmap.begin(); it != urlmap.end(); it++){
		std::cout << "Map key " << (*it).first << " to " << (*it).second << "\n" << std::flush;
	}
	checkAfterRepair = checkAgainAfterRepair;

	//size_t blkid;
	//size_t strpid;
	chunksRepaired = 0;
	currentBlockChecked = 0;
	//finishedRepair = false;
	std::unique_lock<std::mutex> lk(finishedRepairMutex);
	CheckBlock();
	std::cout<<"Out of Check Block of main repair method\n"<< std::flush;
	repairVar.wait(lk, [this]{return this->finishedRepair;});
	lk.unlock();
	CloseAllArchives(0);
	std::cout<<"Out of Archive Closing of main repair method\n"<< std::flush;

	std::unique_ptr<ObjCfg> cfg = std::make_unique<ObjCfg>(objcfg);
		for(size_t i = 0; i < cfg->nbchunks; i++){
			if(redirectionMap.find(cfg->GetDataUrl(i)) != redirectionMap.end()){
				cfg->plgr[i] = redirectionMap[cfg->GetDataUrl(i)];
			}
		}
		for(auto it = cfg->plgr.begin(); it != cfg->plgr.end(); it++){
					std::cout << "Host at " << *it << "\n" << std::flush;
				}
		for(auto it = redirectionMap.begin(); it != redirectionMap.end(); it++)
							std::cout << "Redirect from " << it->first << " to " << it->second << "\n" << std::flush;


		std::cout << "\nCHUNKS REPAIRED: " << chunksRepaired << "\n" << std::flush;
		return cfg;

	//< ID of the block from which we will be reading
		//strpid = (i % objcfg.nbchunks); //< ID of the stripe from which we will be reading
		// generate the file name (blknb/strpnb)
		//std::string fn = objcfg.GetFileName(blkid, strpid);
		// if the block/stripe does not exist it means we are reading passed the end of the file
		//auto itr = urlmap.find(fn);
		//if (itr == urlmap.end()) {


}

void RepairTool::CheckBlock() {
	std::unique_lock<std::mutex> lck(blkmtx);
	redirectMapOffset = 0;
	size_t totalBlocks = 0;
	if (objcfg.nomtfile) {
		totalBlocks = filesize * objcfg.nbchunks / objcfg.blksize;
	} else {
		totalBlocks = lstblk + 1;
	}
	size_t blkid = currentBlockChecked;
	std::cout << "Block " << blkid << " out of " << totalBlocks << "\n"
			<< std::flush;
	if (blkid < totalBlocks) {
		if (!block || block->blkid != blkid)
			block = std::make_shared<block_t>(blkid, reader, objcfg);
		auto blk = block;
		if (!error_correction(blk, this)) {
			currentBlockChecked++;
			lck.unlock();
			CheckBlock();
			std::cout<<"Out of Check Block of Check Block\n"<< std::flush;
		}
	}
	else {
		std::unique_lock<std::mutex> lk(finishedRepairMutex);
		finishedRepair = true;
		std::cout << "Success with repairing!" << std::flush;
		repairVar.notify_all();
	}

}

void RepairTool::WriteChunk(std::shared_ptr<block_t> blk, size_t strpid){
	std::cout << "Writing chunk at stripe " << strpid << "\n" << std::flush;
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
		std::cout << "Data archive " << url << " open\n" << std::flush;
		std::string fn = objcfg.GetFileName( blkid, strpid );
		auto it = writeDataarchs[url]->cdmap.find(fn);
		uint64_t offset = 0;
		if(it != writeDataarchs[url]->cdmap.end()){
			std::cout << "Updating data archive " << url << "\n" << std::flush;
			// the file exists, so we overwrite it
			writeDataarchs[url]->WriteInto(fn, offset, objcfg.chunksize, &(blk->stripes[strpid][0]), nullptr, 0);
			// update chksum??
		}
		else{
			std::cout << "Appending to data archive " << url << ": "<< std::flush;
			for(int u = 0; u < (int)objcfg.chunksize; u++){
				std::cout << blk->stripes[strpid][u];
			}
			std::cout << "\n";
			uint32_t chksum = objcfg.digest(0, &(blk->stripes[strpid][0]), objcfg.chunksize);
			// the file doesnt exist at all, so we append it to the archive (likely a completely new archive)
			writeDataarchs[url]->AppendFile(fn, chksum, objcfg.chunksize, &(blk->stripes[strpid][0]), nullptr, 0);
		}
	}
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

void RepairTool::CloseAllArchives(uint16_t timeout){
	//-------------------------------------------------------------------------
	    // First, check the global status, if we are in an error state just
	    // fail the request.
	    //-------------------------------------------------------------------------
	    //XrdCl::XRootDStatus gst = global_status.get();
	    //if( !gst.IsOK() ) return;// ScheduleHandler( handler, gst );

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
	    XrdCl::Pipeline p2 = ReadSize(0);
	    XrdCl::Async(std::move(p2));

	    std::vector<XrdCl::xattr_t> xav{ {"xrdec.filesize", std::to_string(filesize)},
	                                     {"xrdec.strpver", closeTime.c_str()} };

	    for( size_t i = 0; i < size; ++i )
	    {
	    	// TODO: check whether we also need to close archives that were replaced
	      //-----------------------------------------------------------------------
	      // close ZIP archives with data
	      //-----------------------------------------------------------------------
	      if( writeDataarchs[objcfg.GetDataUrl(i)]->IsOpen() )
	      {
	    	  std::cout << "Data archive with url " << objcfg.GetDataUrl(i) << " is open\n" << std::flush;
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

	    //-------------------------------------------------------------------------
	    // If we were instructed not to create the the additional metadata file
	    // do the simplified close
	    //-------------------------------------------------------------------------
	    if( save_metadata.empty() )
	    {
	      XrdCl::Pipeline p = XrdCl::Parallel( closes ).AtLeast( objcfg.nbchunks );
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

//---------------------------------------------------------------------------
// Factory for creating OpenArchiveImpl objects
//---------------------------------------------------------------------------


void RepairTool::OpenInUpdateMode(XrdCl::ResponseHandler *handler,
		uint16_t timeout) {
	const size_t size = objcfg.plgr.size();
	std::vector<XrdCl::Pipeline> opens;
	opens.reserve(size);
	for (size_t i = 0; i < size; ++i) {
		// generate the URL
		std::string url = objcfg.GetDataUrl(i);
		auto archive = std::make_shared<XrdCl::ZipArchive>(
				Config::Instance().enable_plugins);
		// create the file object
		readDataarchs.emplace(url, archive);
		writeDataarchs.emplace(url, archive);
		// open the archive
		if (objcfg.nomtfile) {
			opens.emplace_back(
					XrdCl::OpenArchive(*readDataarchs[url], url,
							XrdCl::OpenFlags::Update));
		} else
			opens.emplace_back(OpenOnly(*readDataarchs[url], url, true));
	}

	auto pipehndl = [=](const XrdCl::XRootDStatus &st) { // set the central directories in ZIP archives (if we use metadata files)
		auto itr = readDataarchs.begin();
		for (; itr != readDataarchs.end(); ++itr) {
			const std::string &url = itr->first;
			auto &zipptr = itr->second;
			if (zipptr->openstage == XrdCl::ZipArchive::NotParsed)
				zipptr->SetCD(metadata[url]);
			else if (zipptr->openstage != XrdCl::ZipArchive::Done
					&& !metadata.empty())
			{
				// host not reachable or file missing: open completely new archive
				AddMissing(metadata[url], url);
			}
			auto itr = zipptr->cdmap.begin();
			for (; itr != zipptr->cdmap.end(); ++itr) {
				urlmap.emplace(itr->first, url);
				size_t blknb = fntoblk(itr->first);
				if (blknb > lstblk)
					lstblk = blknb;
			}
		}
		metadata.clear();
		// call user handler
		if (handler)
			handler->HandleResponse(new XrdCl::XRootDStatus(st), nullptr);
	};
	// in parallel open the data files and read the metadata
	XrdCl::Pipeline p =
			objcfg.nomtfile ?
					XrdCl::Parallel(opens).AtLeast(objcfg.nbdata) | ReadSize(0)
							| XrdCl::Final(pipehndl) :
					XrdCl::Parallel(ReadMetadata(0),
							XrdCl::Parallel(opens).AtLeast(objcfg.nbdata))
							>> pipehndl;
	XrdCl::Async(std::move(p), timeout);
}

//-----------------------------------------------------------------------
// Add all the entries from given Central Directory to missing and create new archive on other host if necessary
//-----------------------------------------------------------------------
void RepairTool::AddMissing(const buffer_t &cdbuff, const std::string &url) {
	if (objcfg.plgrReplace.size() > 0) {
		const std::string newUrl = objcfg.GetReplacementUrl(0);
		const std::string replacementPlgr = objcfg.plgrReplace[0];
		objcfg.plgrReplace.erase(objcfg.plgrReplace.begin());
		auto newArch = std::make_shared<XrdCl::ZipArchive>(
				Config::Instance().enable_plugins);
		newArch->OpenArchive(newUrl,
				XrdCl::OpenFlags::New | XrdCl::OpenFlags::Write, nullptr, 0);
		writeDataarchs[url] = newArch;
		std::cout << "Creating new archive pointing from url " << url << " to " << newUrl << "\n" << std::flush;
		// save a mapping from old to new urls so user knows where actual data is
		redirectionMap[url] = replacementPlgr;
	}else{
		// TODO: throw error if there's no new url
	}


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
				filesize = std::stoull(size);
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

void RepairTool::Read( size_t blknb, size_t strpnb, buffer_t &buffer, callback_t cb, uint16_t timeout )
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

  std::cout << "Reading from " << url << "\n" << std::flush;

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
                  [zipptr, fn, cb, this]( XrdCl::XRootDStatus &st, XrdCl::ChunkInfo &ch )
                  {
                    //---------------------------------------------------
                    // If read failed there's nothing to do, just pass the
                    // status to user callback
                    //---------------------------------------------------
                    if( !st.IsOK() )
                    {
                      cb( st, 0 );
                      return;
                    }
                    //---------------------------------------------------
                    // Get the checksum for the read data
                    //---------------------------------------------------
                    uint32_t orgcksum = 0;
                    auto s = zipptr->GetCRC32( fn, orgcksum );
                    //---------------------------------------------------
                    // If we cannot extract the checksum assume the data
                    // are corrupted
                    //---------------------------------------------------
                    if( !st.IsOK() )
                    {
                      cb( st, 0 );
                      return;
                    }
                    //---------------------------------------------------
                    // Verify data integrity
                    //---------------------------------------------------
                    uint32_t cksum = objcfg.digest( 0, ch.buffer, ch.length );
                    if( orgcksum != cksum )
                    {
                      cb( XrdCl::XRootDStatus( XrdCl::stError, XrdCl::errDataError ), 0 );
                      return;
                    }
                    //---------------------------------------------------
                    // All is good, we can call now the user callback
                    //---------------------------------------------------
                    cb( XrdCl::XRootDStatus(), ch.length );
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
  if( objcfg.nomtfile && fntoblk( fn ) <= lstblk ) return true;
  // otherwise return false
  return false;
}


}
