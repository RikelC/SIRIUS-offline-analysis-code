#ifndef siriusDataFrame_h
#define siriusDataFrame_h 1
#include "global.hh"
class siriusDataFrame
{
 private:
  
  int metatype ;
  int blocksize ;
  bool endiness ;
  bool blobness ;
  int framesize ;
  uint16_t datasource;
  uint16_t frametype;
  uint8_t rivision;
  uint32_t headersize;
  int  dumpsize;
  int buffersize;

  
  uint32_t eventnumber;
  uint64_t timestamp;
  uint16_t ftrace[TRACE_SIZE];
  int fboardID;
  int fchannelID;
 public:
  siriusDataFrame();
  ~siriusDataFrame();

  void set_framesize(uint32_t i){framesize = i;}
  void set_metatype(int i){metatype = i;}
  void set_blocksize(int i){blocksize = i;}
  void set_blobness(bool i){blobness = i;}
  void set_endiness(bool i){endiness = i;}
  void set_datasource(uint8_t i){datasource = i;}
  void set_frametype(uint16_t i){frametype = i;}
  void set_rivision(uint8_t i){rivision = i;}
  void set_headersize(uint32_t i){headersize =i;}
  void set_dumpsize(int i){dumpsize = i;}
  void set_buffersize(int i){buffersize = i;}
  void set_eventnumber(uint32_t i){eventnumber =i;}
  void set_timestamp(uint64_t i){timestamp =i;}
  void set_boardID(int i){fboardID =i;}
  void set_channelID(int i){fchannelID = i;}
  void set_trace_value(int i, uint16_t v){ftrace[i] = v;}
  //
  uint32_t get_framesize()const {return framesize;}
  int get_metatype()const {return metatype;}
  int get_blocksize()const {return blocksize;}
  bool get_blobness()const {return blobness;}
  bool get_endiness()const {return endiness;}
  uint16_t get_datasource()const {return datasource;}
  uint16_t get_frametype()const {return frametype;}
  uint8_t get_rivision()const {return rivision;}
  uint32_t get_headersize()const {return headersize;}
  int get_dumpsize()const {return dumpsize;}
  int get_bufferzie() const{return buffersize;}
  uint32_t get_eventnumber()const {return eventnumber;}
  uint64_t get_timestamp()const {return timestamp;}
  int get_boardID()const {return fboardID;}
  int get_channelID()const {return fchannelID;}
  uint16_t get_trace_value(int i)const {return ftrace[i];}
  uint16_t* get_trace(){return ftrace;}
  
};
#endif
