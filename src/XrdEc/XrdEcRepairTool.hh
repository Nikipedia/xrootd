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
#ifndef SRC_XRDEC_XRDECREPAIRTOOL_HH_
#define SRC_XRDEC_XRDECREPAIRTOOL_HH_

#include "XrdEc/XrdEcWrtBuff.hh"
#include "XrdEc/XrdEcThreadPool.hh"

#include "XrdEc/XrdEcReader.hh"

#include "XrdCl/XrdClFileOperations.hh"
#include "XrdCl/XrdClParallelOperation.hh"
#include "XrdCl/XrdClZipOperations.hh"

#include "XrdEc/XrdEcObjCfg.hh"

#include "XrdCl/XrdClZipArchive.hh"
#include "XrdCl/XrdClOperations.hh"

#include <string>
#include <unordered_map>
#include <unordered_set>

#include <sys/stat.h>

namespace XrdEc {

struct ThreadEndSemaphore{
	ThreadEndSemaphore(std::shared_ptr<XrdSysSemaphore> s) : sem(s) {}
	~ThreadEndSemaphore(){std::cout << "Semaphore destroyed\n"<<std::flush;sem->Post();sem->Post();}
	std::shared_ptr<XrdSysSemaphore> sem;
};

class RepairTool {

public:
	RepairTool(ObjCfg &objcfg) :
			objcfg(objcfg), reader(objcfg), lstblk(0), filesize(0),
			checkAfterRepair(false){
		currentBlockChecked = 0;
		redirectMapOffset = 0;
		chunksRepaired = 0;
		finishedRepair = false;
	}
	virtual ~RepairTool() {
	}
	void RepairFile(bool checkAgainAfterRepair, XrdCl::ResponseHandler *handler);
	size_t currentBlockChecked;
	uint64_t chunksRepaired;
private:
	void CheckBlock();
	void CheckAllMetadata(std::shared_ptr<ThreadEndSemaphore> sem);
	void InvalidateReplaceArchive(const std::string &url, std::shared_ptr<XrdCl::ZipArchive> zipptr);
	void CompareLFHToCDFH(std::shared_ptr<ThreadEndSemaphore> sem, uint16_t blkid, uint16_t strpid);
	static bool error_correction( std::shared_ptr<block_t> &self, RepairTool *writer );
	static callback_t update_callback(std::shared_ptr<block_t> &self, RepairTool *tool, size_t strpid);
	void OpenInUpdateMode(XrdCl::ResponseHandler *handler,
			uint16_t timeout = 0);
	void CloseAllArchives(XrdCl::ResponseHandler *handler, uint16_t timeout = 0);
	XrdZip::buffer_t GetMetadataBuffer();
	void ReplaceURL(const std::string &url);
	void AddMissing(const buffer_t &cdbuff);
	XrdCl::Pipeline ReadMetadata( size_t index );
    //-----------------------------------------------------------------------
    //! Read size from xattr
    //!
    //! @param index : placement's index
    //-----------------------------------------------------------------------
    XrdCl::Pipeline ReadSize( size_t index );

    std::vector<XrdCl::Pipeline> ReadHealth();

    void Read( size_t blknb, size_t strpnb, buffer_t &buffer, callback_t cb, uint16_t timeout = 0, bool exactControl = false);

    bool IsMissing(const std::string &fn);

    //-----------------------------------------------------------------------
    //! Parse metadata from chunk info object
    //!
    //! @param ch : chunk info object returned by a read operation
    //-----------------------------------------------------------------------
    bool ParseMetadata( XrdCl::ChunkInfo &ch );
    void WriteChunk(std::shared_ptr<block_t> blk, size_t strpid);

	ObjCfg &objcfg;
	// unused reader only for initialization of block_t
	Reader reader;
	std::vector<std::shared_ptr<XrdCl::File>>        metadataarchs;      //< ZIP archives with metadata
	Reader::dataarchs_t readDataarchs; //> map URL to ZipArchive object
	Reader::dataarchs_t writeDataarchs; //> map URL to ZipArchives for writing (may be different!)
	Reader::metadata_t metadata;  //> map URL to CD metadata
	Reader::urlmap_t urlmap;    //> map blknb/strpnb (data chunk) to URL
	Reader::urlmap_t redirectionMap;
	Reader::missing_t missing;   //> set of missing stripes
	std::shared_ptr<block_t> block;  //> cache for the block we are reading from
	std::mutex blkmtx;    //> mutex guarding the block from parallel access
	size_t lstblk;    //> last block number
	uint64_t filesize;  //> file size (obtained from xattr)

	// for replacing URLs
	std::mutex urlMutex;


    int redirectMapOffset;

	std::mutex finishedRepairMutex;
	std::condition_variable repairVar;
	bool finishedRepair;

	bool checkAfterRepair;
};
}

#endif
