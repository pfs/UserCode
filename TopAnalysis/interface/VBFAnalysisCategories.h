#ifndef vbfanalysiscategories_h
#define vbfanalysiscategories_h
#include <iostream>
using namespace std;
namespace vbf {
  TString categoryNames[]={"LowVPtLowMJJ","LowVPtHighMJJ","HighVPtLowMJJ","HighVPtHighMJJ","HighVPt","LowVPt"," LooseVPt"};
  struct Category{
    float EE,MM,A,LowVPt,HighVPt,LowMJJ,HighMJJ,AllMJJ, LooseVPt;
    Category(){ reset(); }
    Category(std::vector<bool> &cat){
      reset();
      set(cat);
    };  
    void reset(){
      std::vector<bool> cat(8,false);
      set(cat);
    };
    void set(std::vector<bool> &cat){
      EE              = (float) cat[0];
      MM              = (float) cat[1];
      A               = (float) cat[2];
      LowVPt          = (float) cat[3];
      HighVPt         = (float) cat[4];
      LowMJJ          = (float) cat[5];
      HighMJJ         = (float) cat[6];
      AllMJJ          = (float) cat[7];
      LooseVPt        = (float) cat[8];

    };
    std::vector<TString> getChannelTags() {
      std::vector<TString> chTags;
      TString chTag("");    
      if(EE>0) chTag="EE";
      if(MM>0) chTag="MM";
      if(A>0)  chTag="A";
      if(chTag=="") return chTags;
      chTags.push_back(chTag);
      if(LowVPt>0  && LowMJJ>0)                  {chTags.push_back(categoryNames[0]+chTag); }//cout<<categoryNames[0]<<endl;}
      if(LowVPt>0  && HighMJJ>0)                 {chTags.push_back(categoryNames[1]+chTag); }//cout<<categoryNames[1]<<endl;}
      if(HighVPt>0 && LowMJJ>0)                  {chTags.push_back(categoryNames[2]+chTag); }//cout<<categoryNames[2]<<endl;}
      if(HighVPt>0 && HighMJJ>0)                 {chTags.push_back(categoryNames[3]+chTag); }//cout<<categoryNames[3]<<endl;}
      if(HighVPt>0 && AllMJJ>0)                  {chTags.push_back(categoryNames[4]+chTag); }//cout<<categoryNames[4]<<endl;}
      if(LowVPt>0 && AllMJJ>0)                   {chTags.push_back(categoryNames[5]+chTag); }//cout<<categoryNames[5]<<endl;}
      if(LowVPt>0 && LooseVPt>0)                 {chTags.push_back(categoryNames[6]+chTag); }//cout<<categoryNames[5]<<endl;}
      if(HighVPt>0 && LooseVPt>0)                 {chTags.push_back(categoryNames[7]+chTag); }//cout<<categoryNames[5]<<endl;} 
      return chTags;
    }
  };

}

#endif
