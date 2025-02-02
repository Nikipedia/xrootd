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

#include <cppunit/extensions/HelperMacros.h>
#include "TestEnv.hh"
#include "CppUnitXrdHelpers.hh"

#include "XrdEc/XrdEcStrmWriter.hh"
#include "XrdEc/XrdEcReader.hh"
#include "XrdEc/XrdEcObjCfg.hh"
#include "XrdEc/XrdEcRepairTool.hh"

#include "XrdCl/XrdClMessageUtils.hh"

#include "XrdZip/XrdZipCDFH.hh"

#include <string>
#include <memory>
#include <limits>

#include <unistd.h>
#include <cstdio>
#include <ftw.h>
#include <sys/stat.h>
#include <sys/types.h>

using namespace XrdEc;

//------------------------------------------------------------------------------
// Declaration
//------------------------------------------------------------------------------
class MicroTest: public CppUnit::TestCase
{
  public:
    CPPUNIT_TEST_SUITE( MicroTest );
      CPPUNIT_TEST( AlignedWriteTest );
      CPPUNIT_TEST( SmallWriteTest );
      CPPUNIT_TEST( BigWriteTest );
      CPPUNIT_TEST( AlignedWrite1MissingTest );
      CPPUNIT_TEST( AlignedWrite2MissingTest );
      CPPUNIT_TEST( AlignedWriteTestIsalCrcNoMt );
      CPPUNIT_TEST( SmallWriteTestIsalCrcNoMt );
      CPPUNIT_TEST( BigWriteTestIsalCrcNoMt );
      CPPUNIT_TEST( AlignedWrite1MissingTestIsalCrcNoMt );
      CPPUNIT_TEST( AlignedWrite2MissingTestIsalCrcNoMt );
      CPPUNIT_TEST( AlignedRepairNoHostTest );
      CPPUNIT_TEST( AlignedRepairOneChunkTest);
      CPPUNIT_TEST( RandomizedRepairTest );
      CPPUNIT_TEST( RepairNoCorruptionTest );
      CPPUNIT_TEST( RandomizedDoubleRepairTest );
    CPPUNIT_TEST_SUITE_END();

    int testedChunkCount;
	int failedReads;

    void Init( bool usecrc32c );

    void InitRepair(bool usecrc32c);

    inline void AlignedWriteTestImpl( bool usecrc32c )
    {
      // create the data and stripe directories
      Init( usecrc32c );
      // run the test
      AlignedWriteRaw();
      // verify that we wrote the data correctly
      Verify(true);
      // clean up the data directory
      CleanUp();
    }

    inline void AlignedWriteTest()
    {
      AlignedWriteTestImpl( true );
    }

    inline void AlignedWriteTestIsalCrcNoMt()
    {
      AlignedWriteTestImpl( false );
    }

    /*
     * corruption type: 0 = missing host, 1 = single corrupt, 2 = random corrupt
     */
    inline void AlignedRepairTestImpl(bool usecrc32c, int corruptionType, bool mustHaveErrors = true){
    	std::cout<<"Repair Test started\n\n\n\n"<<std::flush;
    	uint64_t seed = std::chrono::system_clock::now().time_since_epoch().count();
    	//last error:
    	//uint64_t seed = 1659356007585588302;
    	randomSeed = seed;

    	InitRepair(usecrc32c);
    	uint64_t seed2 = std::chrono::system_clock::now().time_since_epoch().count();
    	//last error:
    	//uint64_t seed2= 1659356007585918572;

    	std::cout << "Random Seed FileGen: " << seed << " Random Seed Corruption: " << seed2<<"\n" << std::flush;

		// run the test

    	AlignedRandomWriteRaw(seed);

    	Verify(false);

    	switch(corruptionType){
    	case 0: UrlNotReachable(4); break;
    	case 1: CorruptChunk(1,1); break;
    	case 2: CorruptRandom(seed2, 1); break;
    	case 3: CorruptRandom(seed2, 2); break;
    	default: break;
    	}
    	if(mustHaveErrors)
    		VerifyAnyErrorExists();

    	auto oldPlgr = std::vector<std::string>(objcfg->plgr);

		XrdEc::RepairTool r(*objcfg);
		// shouldnt affect the repair tool's objcfg object since it was copied
		//objcfg.reset();
		//std::shared_ptr<ObjCfg> updatedConfig = r.RepairFile(false, nullptr);
		r.RepairFile(false, nullptr);
		//objcfg = std::make_unique<ObjCfg>(*updatedConfig);
		// shouldn't affect the updated objcfg since its a copy
		//updatedConfig.reset();
		if(mustHaveErrors)
			CPPUNIT_ASSERT(r.chunksRepaired > 0);

		std::cout<<"Repaired file, starting verification\n"<<std::flush;

		Verify(false);
		// clean up
		if(corruptionType == 0){
			objcfg->plgr = oldPlgr;
			UrlReachable(4);
		//UrlReachable(3);
		//UrlReachable(4);
		}
		CleanUp();
    }

    inline void AlignedRepairNoHostTest(){
    	AlignedRepairTestImpl(true, 0);
    }

    inline void AlignedRepairOneChunkTest(){
    	AlignedRepairTestImpl(true, 1);
    }

    inline void RandomizedRepairTest(){
    	AlignedRepairTestImpl(true, 2);
    }

    inline void RepairNoCorruptionTest(){
    	AlignedRepairTestImpl(true, -1, false);
    }

    inline void RandomizedDoubleRepairTest(){
    	AlignedRepairTestImpl(true, 3);
    }

    inline void AlignedWrite1MissingTestImpl( bool usecrc32c )
    {
      // initialize directories
      Init( usecrc32c );
      UrlNotReachable( 2 );
      // run the test
      AlignedWriteRaw();
      // verify that we wrote the data correctly
      Verify(true);
      // clean up
      UrlReachable( 2 );
      CleanUp();
    }

    inline void AlignedWrite1MissingTest()
    {
      AlignedWrite1MissingTestImpl( true );
    }

    inline void AlignedWrite1MissingTestIsalCrcNoMt()
    {
      AlignedWrite1MissingTestImpl( false );
    }

    inline void AlignedWrite2MissingTestImpl( bool usecrc32c )
    {
      // initialize directories
      Init( usecrc32c );
      UrlNotReachable( 2 );
      UrlNotReachable( 3 );
      // run the test
      AlignedWriteRaw();
      // verify that we wrote the data correctly
      Verify(true);
      // clean up
      UrlReachable( 2 );
      UrlReachable( 3 );
      CleanUp();
    }

    inline void AlignedWrite2MissingTest()
    {
      AlignedWrite2MissingTestImpl( true );
    }

    inline void AlignedWrite2MissingTestIsalCrcNoMt()
    {
      AlignedWrite2MissingTestImpl( false );
    }

    void VarlenWriteTest( uint32_t wrtlen, bool usecrc32c );

    inline void SmallWriteTest()
    {
      VarlenWriteTest( 7, true );
    }

    inline void SmallWriteTestIsalCrcNoMt()
    {
      VarlenWriteTest( 7, false );
    }

    void BigWriteTest()
    {
      VarlenWriteTest( 77, true );
    }

    void BigWriteTestIsalCrcNoMt()
    {
      VarlenWriteTest( 77, false );
    }

    void Verify(bool repairAllow)
    {
      ReadVerifyAll(repairAllow);
      if(repairAllow)
    	  CorruptedReadVerify();
    }

    void VerifyAnyErrorExists();

    void CleanUp();

    inline void ReadVerifyAll(bool repairAllow)
    {
      AlignedReadVerify(repairAllow);
      std::cout<<"Aligned Read success\n"<<std::flush;
      PastEndReadVerify(repairAllow);
      std::cout<<"Past End Read success\n"<<std::flush;
      SmallChunkReadVerify(repairAllow);
      std::cout<<"Small Read success\n"<<std::flush;
      BigChunkReadVerify(repairAllow);
      std::cout<<"Big Read success\n"<<std::flush;

      for( size_t i = 0; i < 10; ++i )
        RandomReadVerify(repairAllow);
      std::cout << "All reads verified"<<std::flush;
    }

    void ReadVerify( uint32_t rdsize, bool repairAllow, uint64_t maxrd = std::numeric_limits<uint64_t>::max() );
    void RandomReadVerify(bool repairAllow);

    void Corrupted1stBlkReadVerify(bool repairAllow);

    inline void AlignedReadVerify(bool repairAllow)
    {
      ReadVerify( chsize, repairAllow, rawdata.size() );
    }

    inline void PastEndReadVerify(bool repairAllow)
    {
      ReadVerify( chsize, repairAllow );
    }

    inline void SmallChunkReadVerify(bool repairAllow)
    {
      ReadVerify( 5, repairAllow );
    }

    inline void BigChunkReadVerify(bool repairAllow)
    {
      ReadVerify( 23 , repairAllow);
    }

    /*
     * switches off some hosts, so this test is only useful with error correction allowed
     */
    void CorruptedReadVerify();

    void CorruptChunk( size_t blknb, size_t strpnb );

    void CorruptRandom(uint64_t seed, uint32_t numCorruptions);

    void UrlNotReachable( size_t index );
    void UrlReachable( size_t index );




  private:

    void AlignedWriteRaw();

    void AlignedRandomWriteRaw(uint64_t seed);

    static callback_t read_callback(MicroTest *self, size_t blkid, size_t strpid);

    void copy_rawdata( char *buffer, size_t size )
    {
      const char *begin = buffer;
      const char *end   = begin + size;
      std::copy( begin, end, std::back_inserter( rawdata ) );
    }



    std::string datadir;
    std::unique_ptr<ObjCfg> objcfg;
    std::vector<std::string> replaceHosts = {"host1", "host2", "host3"};

    static const size_t nbdata   = 4;
    static const size_t nbparity = 2;
    //static const size_t chsize   = 16;
    static const size_t chsize   = 128;
    static const size_t nbiters  = 16;

    static const size_t lfhsize  = 30;

    uint32_t randomSeed;

    std::vector<char> rawdata;
};

CPPUNIT_TEST_SUITE_REGISTRATION( MicroTest );


void MicroTest::Init( bool usecrc32c )
{
  objcfg.reset( new ObjCfg( "test.txt", nbdata, nbparity, chsize, usecrc32c, true ) );
  rawdata.clear();

  char cwdbuff[1024];
  char *cwdptr = getcwd( cwdbuff, sizeof( cwdbuff ) );
  CPPUNIT_ASSERT( cwdptr );
  std::string cwd = cwdptr;
  // create the data directory
  datadir = cwd + "/data";
  CPPUNIT_ASSERT( mkdir( datadir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH ) == 0 );
  // create a directory for each stripe
  size_t nbstrps = objcfg->nbdata + 2 * objcfg->nbparity;
  for( size_t i = 0; i < nbstrps; ++i )
  {
    std::stringstream ss;
    ss << std::setfill('0') << std::setw( 2 ) << i;
    std::string strp = datadir + '/' + ss.str() + '/';
    objcfg->plgr.emplace_back( strp );
    CPPUNIT_ASSERT( mkdir( strp.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH ) == 0 );
  }
}

void MicroTest::InitRepair( bool usecrc32c )
{
  objcfg.reset( new ObjCfg( "test.txt", nbdata, nbparity, chsize, usecrc32c, true ) );
  rawdata.clear();

  char cwdbuff[1024];
  char *cwdptr = getcwd( cwdbuff, sizeof( cwdbuff ) );
  CPPUNIT_ASSERT( cwdptr );
  std::string cwd = cwdptr;
  // create the data directory
  datadir = cwd + "/data";
  CPPUNIT_ASSERT( mkdir( datadir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH ) == 0 );
  // create a directory for each stripe
  //size_t nbstrps = objcfg->nbdata + 2 * objcfg->nbparity;
  size_t nbstrps = objcfg->nbdata + objcfg->nbparity;
  for( size_t i = 0; i < nbstrps; ++i )
  {
    std::stringstream ss;
    ss << std::setfill('0') << std::setw( 2 ) << i;
    std::string strp = datadir + '/' + ss.str() + '/';
    objcfg->plgr.emplace_back( strp );
    CPPUNIT_ASSERT( mkdir( strp.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH ) == 0 );
  }
  for(size_t i = 0; i < replaceHosts.size(); i++){
	    std::string strp = datadir + '/' + replaceHosts[i] + '/';
	  objcfg->plgrReplace.emplace_back(strp);
	  CPPUNIT_ASSERT( mkdir( strp.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH ) == 0 );
  }
}


void MicroTest::CorruptChunk( size_t blknb, size_t strpnb )
{
  Reader reader( *objcfg );
  // open the data object
  XrdCl::SyncResponseHandler handler1;
  reader.Open( &handler1 );
  handler1.WaitForResponse();
  XrdCl::XRootDStatus *status = handler1.GetStatus();
  CPPUNIT_ASSERT_XRDST( *status );
  delete status;

  // get the CD buffer
  std::string fn     = objcfg->GetFileName( blknb, strpnb );
  std::string url    = reader.urlmap[fn];
  buffer_t    cdbuff = reader.dataarchs[url]->GetCD();

  // close the data object
  XrdCl::SyncResponseHandler handler2;
  reader.Close( &handler2 );
  handler2.WaitForResponse();
  status = handler2.GetStatus();
  CPPUNIT_ASSERT_XRDST( *status );
  delete status;

  // parse the CD buffer
  const char *buff   = cdbuff.data();
  size_t      size   = cdbuff.size();
  XrdZip::cdvec_t cdvec;
  XrdZip::cdmap_t cdmap;
  std::tie(cdvec, cdmap ) = XrdZip::CDFH::Parse( buff, size );

  // now corrupt the chunk (put wrong checksum)
  XrdZip::CDFH &cdfh = *cdvec[cdmap[fn]];
  uint64_t offset = cdfh.offset + lfhsize + fn.size(); // offset of the data
  XrdCl::File f;
  XrdCl::XRootDStatus status2 = f.Open( url, XrdCl::OpenFlags::Write );
  CPPUNIT_ASSERT_XRDST( status2 );
  std::string str = "XXXXXXXX";
  status2 = f.Write( offset, str.size(), str.c_str() );
  CPPUNIT_ASSERT_XRDST( status2 );
  status2 = f.Close();
  CPPUNIT_ASSERT_XRDST( status2 );
}

void MicroTest::CorruptRandom(uint64_t seed, uint32_t numCorruptions){
	static std::default_random_engine random_engine(seed);
	std::uniform_int_distribution<uint32_t> hostToCorrupt(0, nbdata+nbparity-1);
	std::vector<uint64_t> hosts;
	std::vector<uint64_t> maxCorruptSizes;
	for(size_t i = 0; i < numCorruptions; i++){
		uint64_t host = hostToCorrupt(random_engine);
		hosts.emplace_back(host);
	}
	Reader reader(*objcfg);
	// open the data object
	XrdCl::SyncResponseHandler handler1;
	reader.Open(&handler1);
	handler1.WaitForResponse();
	XrdCl::XRootDStatus *status = handler1.GetStatus();
	CPPUNIT_ASSERT_XRDST(*status);
	delete status;
	for(size_t h = 0; h < nbdata+nbparity; h++){
			// get the CD buffer
			std::string url = objcfg->GetDataUrl(h);
			auto zipptr = reader.dataarchs[url];
			uint64_t cdOffset =
					zipptr->zip64eocd ?
							zipptr->zip64eocd->cdOffset : zipptr->eocd->cdOffset;
			uint64_t cdLength = zipptr->zip64eocd? zipptr->zip64eocd->cdSize:zipptr->eocd->cdSize;
			uint64_t eocdLength = XrdZip::EOCD::eocdBaseSize;

			maxCorruptSizes.emplace_back(cdOffset + cdLength + eocdLength);
	}
	// close the data object
	XrdCl::SyncResponseHandler handler2;
	reader.Close(&handler2);
	handler2.WaitForResponse();
	status = handler2.GetStatus();
	CPPUNIT_ASSERT_XRDST(*status);
	delete status;

	for(size_t i = 0; i < numCorruptions; i++){
		std::string url = objcfg->GetDataUrl(hosts[i]);
		std::stringstream ss;
		ss << "cp " << url << " " << url << "_noncorrupt";
		std::string s = ss.str();
		system(s.data());

	// now corrupt the host
	XrdCl::File f;
	XrdCl::XRootDStatus status2 = f.Open(url, XrdCl::OpenFlags::Update);
	CPPUNIT_ASSERT_XRDST(status2);


	std::uniform_int_distribution<uint32_t> lengthRandom(1, chsize);
	uint64_t size = lengthRandom(random_engine);
	uint64_t writeSize = size;

	std::uniform_int_distribution<uint32_t> offsetRandom(0, maxCorruptSizes[hosts[i]] - size);
	uint64_t offset = offsetRandom(random_engine);


	std::uniform_int_distribution<uint32_t> letter( 0, 25 );

	std::vector<char> buffer;

	std::cout << "Corrupting host " << url << " at offset " << offset << " with size " << size << ":\n" << std::flush;

	while(size > 0){
		uint32_t letterAdd = letter(random_engine);
		buffer.push_back('A' + letterAdd);
		std::cout << (char)('A'+letterAdd);
		size -= 1;
	}
	std::cout << "\n" << std::flush;
	status2 = f.Write(offset, writeSize, buffer.data());
	CPPUNIT_ASSERT_XRDST(status2);
	status2 = f.Close();
	CPPUNIT_ASSERT_XRDST(status2);

	}
	std::cout << "Successfully corrupted archive\n" << std::flush;
}

void MicroTest::UrlNotReachable( size_t index )
{
  XrdCl::URL url( objcfg->plgr[index] );
  CPPUNIT_ASSERT( chmod( url.GetPath().c_str(), 0 ) == 0 );
}

void MicroTest::UrlReachable( size_t index )
{
  XrdCl::URL url( objcfg->plgr[index] );
  mode_t mode = S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH |
                S_IXUSR | S_IXGRP | S_IXOTH;
  CPPUNIT_ASSERT( chmod( url.GetPath().c_str(), mode ) == 0 );
}

void MicroTest::CorruptedReadVerify()
{
  UrlNotReachable( 0 );
  ReadVerifyAll(true);
  UrlNotReachable( 1 );
  ReadVerifyAll(true);
  UrlReachable( 0 );
  UrlReachable( 1 );

  CorruptChunk( 0, 0 );
  ReadVerifyAll(true);

  CorruptChunk( 0, 1 );
  ReadVerifyAll(true);

  CorruptChunk( 0, 2 );
  Corrupted1stBlkReadVerify(true);
}

callback_t MicroTest::read_callback(MicroTest *self, size_t blkid, size_t strpid) {
	return [self, blkid, strpid](const XrdCl::XRootDStatus &st,
			uint32_t) mutable {
			if (st.IsOK()){
				self->testedChunkCount++;
				//std::cout << "successful read of block " << (int)blkid << ", Stripe "<<(int)strpid<< "\n"<<std::flush;
			}
			else{
				self->failedReads++;
				std::cout << "Fail in block " << (int)blkid << ", Stripe "<<(int)strpid<<": "<< st.code << ": "<< st.GetErrorMessage()<< "\n" << std::flush;
				if(st.code == XrdCl::errCorruptedHeader){
					std::cout << "Corrupted metadata in block " << (int)blkid << ", Stripe "<<(int)strpid<< "\n" << std::flush;
				}
			}
		};
}

void MicroTest::VerifyAnyErrorExists(){
	RepairTool tool(*objcfg);
	XrdCl::SyncResponseHandler handlerRepair;
	tool.CheckFile(&handlerRepair);
	handlerRepair.WaitForResponse();
	XrdCl::XRootDStatus *status = handlerRepair.GetStatus();
		if(!status->IsOK()){
			std::cout << "Archives can't be opened as desired";
			return;
		}

	Reader reader(*objcfg);
	// open the data object
	XrdCl::SyncResponseHandler handler1;
	reader.Open(&handler1);
	handler1.WaitForResponse();
	status = handler1.GetStatus();
	//CPPUNIT_ASSERT_XRDST(*status);
	if(!status->IsOK()){
		std::cout << "Archives can't be opened as desired";
		return;
	}

	std::cout << "Opened the reader successfully\n"<< std::flush;
	uint64_t rdsize = chsize;
	uint64_t maxrd = rawdata.size();
	uint64_t rdoff = 0;
	char *rdbuff = new char[rdsize];
	uint32_t bytesrd = 0;
	uint64_t total_bytesrd = 0;
	bool errorFound = false;
	//std::cout << "Reading with offset ";
	do {
		XrdCl::SyncResponseHandler h;
		//std::cout << rdoff << std::flush;
		reader.Read(rdoff, rdsize, rdbuff, &h, 0);
		h.WaitForResponse();
		status = h.GetStatus();
		if(!status->IsOK()){
			std::cout << "Couldn't read " << (int)rdsize << " at offset " << (int) rdoff << "\n"<< std::flush;
			errorFound = true;
			delete status;
			break;
		}
		// get the actual result
		auto rsp = h.GetResponse();
		XrdCl::ChunkInfo *ch = nullptr;
		rsp->Get(ch);
		bytesrd = ch->length;
		std::string result(reinterpret_cast<char*>(ch->buffer), bytesrd);
		// get the expected result
		size_t rawoff = rdoff;
		size_t rawsz = rdsize;
		if (rawoff + rawsz > rawdata.size())
			rawsz = rawdata.size() - rawoff;
		std::string expected(rawdata.data() + rawoff, rawsz);
		std::cout << result << " - " << expected << "\n" << std::flush;

		// make sure the expected and actual results are the same
		if(result != expected){
			errorFound = true;
		}
		delete status;
		delete rsp;
		rdoff += bytesrd;
		total_bytesrd += bytesrd;
	} while (bytesrd == rdsize && total_bytesrd < maxrd && !errorFound);
	delete[] rdbuff;

	// close the data object
	XrdCl::SyncResponseHandler handler2;
	reader.Close(&handler2);
	handler2.WaitForResponse();
	status = handler2.GetStatus();
	CPPUNIT_ASSERT_XRDST(*status);



	delete status;

	CPPUNIT_ASSERT(errorFound);
}



void MicroTest::ReadVerify( uint32_t rdsize, bool repairAllow, uint64_t maxrd )
{
	std::cout << "Read Verification started\n" << std::flush;

	// check metadata and checksums
	RepairTool tool(*objcfg);
	XrdCl::SyncResponseHandler handlerRepair;
	tool.CheckFile(&handlerRepair);
	handlerRepair.WaitForResponse();
	XrdCl::XRootDStatus *status = handlerRepair.GetStatus();
	CPPUNIT_ASSERT_XRDST( *status );

  Reader reader( *objcfg );
  // open the data object
  XrdCl::SyncResponseHandler handler1;
  reader.Open( &handler1 );
  handler1.WaitForResponse();
  status = handler1.GetStatus();
  CPPUNIT_ASSERT_XRDST( *status );
  delete status;
  
  uint64_t  rdoff  = 0;
  char     *rdbuff = new char[rdsize]; 
  uint32_t  bytesrd = 0;
  uint64_t  total_bytesrd = 0;
  //std::cout << "Reading with offset ";
  do
  {
    XrdCl::SyncResponseHandler h;
    //std::cout << rdoff << std::flush;
    reader.Read( rdoff, rdsize, rdbuff, &h, 0);
    h.WaitForResponse();
    status = h.GetStatus();
    CPPUNIT_ASSERT_XRDST( *status );
    // get the actual result
    auto rsp = h.GetResponse();
    XrdCl::ChunkInfo *ch = nullptr;
    rsp->Get( ch );
    bytesrd = ch->length;
    std::string result( reinterpret_cast<char*>( ch->buffer ), bytesrd );
    // get the expected result
    size_t rawoff = rdoff;
    size_t rawsz  = bytesrd;
    if( rawoff + rawsz > rawdata.size() ) rawsz = rawdata.size() - rawoff;
    std::string expected( rawdata.data() + rawoff, rawsz );
    // make sure the expected and actual results are the same
    if(result != expected){
    	std::cout << "Error at offset "<< (int)rdoff << " data doesn't match: "<<expected << " but got: " << result << "\n" << std::flush;
    	std::cout << "Sizes: Expected has size " << expected.size() << ", result has " << result.size() << "\n" << std::flush;
    }
    CPPUNIT_ASSERT( result == expected );
    delete status;
    delete rsp;
    rdoff += bytesrd;
    total_bytesrd += bytesrd;
  }
  while( bytesrd == rdsize && total_bytesrd < maxrd );
  delete[] rdbuff;
 
  // close the data object
  XrdCl::SyncResponseHandler handler2;
  reader.Close( &handler2 );
  handler2.WaitForResponse();
  status = handler2.GetStatus();
  CPPUNIT_ASSERT_XRDST( *status );
  delete status;
}

void MicroTest::RandomReadVerify(bool repairAllow)
{
  size_t filesize = rawdata.size();
  // better reproducibility: set seed
  static std::default_random_engine random_engine( randomSeed );
  std::uniform_int_distribution<uint32_t> offdistr( 0, filesize );
  uint64_t rdoff = offdistr( random_engine );
  std::uniform_int_distribution<uint32_t> lendistr( rdoff, filesize + 32 );
  uint32_t rdlen = lendistr( random_engine );

  Reader reader( *objcfg );
  // open the data object
  XrdCl::SyncResponseHandler handler1;
  reader.Open( &handler1 );
  handler1.WaitForResponse();
  XrdCl::XRootDStatus *status = handler1.GetStatus();
  CPPUNIT_ASSERT_XRDST( *status );
  delete status;

  // read the data
  char *rdbuff = new char[rdlen];
  XrdCl::SyncResponseHandler h;
  reader.Read( rdoff, rdlen, rdbuff, &h, 0);
  h.WaitForResponse();
  status = h.GetStatus();
  CPPUNIT_ASSERT_XRDST( *status );
  // get the actual result
  auto rsp = h.GetResponse();
  XrdCl::ChunkInfo *ch = nullptr;
  rsp->Get( ch );
  uint32_t bytesrd = ch->length;
  std::string result( reinterpret_cast<char*>( ch->buffer ), bytesrd );
  // get the expected result
  size_t rawoff = rdoff;
  size_t rawlen  = rdlen;
  if( rawoff > rawdata.size() ) rawlen = 0;
  else if( rawoff + rawlen > rawdata.size() ) rawlen = rawdata.size() - rawoff;
  std::string expected( rawdata.data() + rawoff, rawlen );
  // make sure the expected and actual results are the same
  CPPUNIT_ASSERT( result == expected );
  delete status;
  delete rsp;
  delete[] rdbuff;

  // close the data object
  XrdCl::SyncResponseHandler handler2;
  reader.Close( &handler2 );
  handler2.WaitForResponse();
  status = handler2.GetStatus();
  CPPUNIT_ASSERT_XRDST( *status );
  delete status;
}

void MicroTest::Corrupted1stBlkReadVerify(bool repairAllow)
{
  uint64_t rdoff = 0;
  uint32_t rdlen = objcfg->datasize;

  Reader reader( *objcfg );
  // open the data object
  XrdCl::SyncResponseHandler handler1;
  reader.Open( &handler1 );
  handler1.WaitForResponse();
  XrdCl::XRootDStatus *status = handler1.GetStatus();
  CPPUNIT_ASSERT_XRDST( *status );
  delete status;

  // read the data
  char *rdbuff = new char[rdlen];
  XrdCl::SyncResponseHandler h;
  reader.Read( rdoff, rdlen, rdbuff, &h, 0 );
  h.WaitForResponse();
  status = h.GetStatus();
  CPPUNIT_ASSERT( status->status == XrdCl::stError &&
                 status->code == XrdCl::errDataError );
  delete status;
  delete[] rdbuff;

  // close the data object
  XrdCl::SyncResponseHandler handler2;
  reader.Close( &handler2 );
  handler2.WaitForResponse();
  status = handler2.GetStatus();
  CPPUNIT_ASSERT_XRDST( *status );
  delete status;
}

int unlink_cb(const char *fpath, const struct stat *sb, int typeflag, struct FTW *ftwbuf)
{
  int rc = remove( fpath );
  CPPUNIT_ASSERT( rc == 0 );
  return rc;
}

void MicroTest::CleanUp()
{
  // delete the data directory
  nftw( datadir.c_str(), unlink_cb, 64, FTW_DEPTH | FTW_PHYS );
}

void MicroTest::AlignedWriteRaw()
{
  char buffer[objcfg->chunksize];
  StrmWriter writer( *objcfg );
  // open the data object
  XrdCl::SyncResponseHandler handler1;
  writer.Open( &handler1 );
  handler1.WaitForResponse();
  XrdCl::XRootDStatus *status = handler1.GetStatus();
  CPPUNIT_ASSERT_XRDST( *status );
  delete status;
  // write to the data object
  for( size_t i = 0; i < nbiters; ++i )
  {
    memset( buffer, 'A' + i, objcfg->chunksize );
    writer.Write( objcfg->chunksize, buffer, nullptr );
    copy_rawdata( buffer, sizeof( buffer ) );
  }
  XrdCl::SyncResponseHandler handler2;
  writer.Close( &handler2 );
  handler2.WaitForResponse();
  status = handler2.GetStatus();
  CPPUNIT_ASSERT_XRDST( *status );
  delete status;
}

void MicroTest::AlignedRandomWriteRaw(uint64_t seed)
{

  StrmWriter writer( *objcfg );
  // open the data object
  XrdCl::SyncResponseHandler handler1;
  writer.Open( &handler1 );
  handler1.WaitForResponse();
  XrdCl::XRootDStatus *status = handler1.GetStatus();
  CPPUNIT_ASSERT_XRDST( *status );
  delete status;
  // generate random stuff
  static std::default_random_engine random_engine( seed );
  std::uniform_int_distribution<uint32_t> filelength( 256, 2560 );
  // data should be equally distributable among all nodes?
  uint64_t size = filelength( random_engine );
  char buffer[size];
  std::uniform_int_distribution<uint32_t> letter( 0, 25 );

  for( size_t i = 0; i < size; ++i )
  {
	uint64_t randLet = letter( random_engine );
    memset( buffer + i, 'A' + randLet, 1 );

  }
  writer.Write( size, buffer, nullptr );
  copy_rawdata( buffer, sizeof( buffer ) );
  XrdCl::SyncResponseHandler handler2;
  writer.Close( &handler2 );
  handler2.WaitForResponse();
  status = handler2.GetStatus();
  CPPUNIT_ASSERT_XRDST( *status );
  delete status;
}

void MicroTest::VarlenWriteTest( uint32_t wrtlen, bool usecrc32c )
{
  // create the data and stripe directories
  Init( usecrc32c );
  // open the data object
  StrmWriter writer( *objcfg );
  XrdCl::SyncResponseHandler handler1;
  writer.Open( &handler1 );
  handler1.WaitForResponse();
  XrdCl::XRootDStatus *status = handler1.GetStatus();
  CPPUNIT_ASSERT_XRDST( *status );
  delete status;
  // write the data
  char     wrtbuff[wrtlen];
  size_t   bytesleft = nbiters * objcfg->chunksize;
  size_t   i = 0;
  while( bytesleft > 0 )
  {
    if( wrtlen > bytesleft ) wrtlen = bytesleft;
    memset( wrtbuff, 'A' + i, wrtlen );
    writer.Write( wrtlen, wrtbuff, nullptr );
    copy_rawdata( wrtbuff, wrtlen );
    bytesleft -= wrtlen;
    ++i;
  }
  XrdCl::SyncResponseHandler handler2;
  writer.Close( &handler2 );
  handler2.WaitForResponse();
  status = handler2.GetStatus();
  CPPUNIT_ASSERT_XRDST( *status );
  delete status;

  // verify that we wrote the data correctly
  Verify(true);
  // clean up the data directory
  CleanUp();
}

