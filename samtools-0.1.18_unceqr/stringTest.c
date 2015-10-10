#include <stdio.h>
#include <string.h>


typedef char * string;
int tti;
int ttj;


//s= "aa";



void main(void){
  //putchar('a');

  string s[10];
  s[0] = "sss";
  s[1] = "pp";
  s[2] = "pp";
  s[3] = "sss";
  s[4] = "23";
  s[5] = "23";

  string insAllele[100];
  int insAlleleCnt[100];
  int sf=0;
  int flag=0;

 for(tti=0;tti<6;tti++){ // for each of strings
    flag=0;
    for(ttj=0;ttj<sf;ttj++){	
      //if(strcmp(insAllele[ttj],s[tti])>0){
      if( strcmp(insAllele[ttj],s[tti])==0 ){
        insAlleleCnt[ttj]++;
        flag=1;
        break;
      }
    }
    if(flag) continue;
    insAllele[sf]=s[tti];
    insAlleleCnt[sf] = 1;
    sf++;
  }



  //print it out
  for(tti=0;tti<sf;tti++){
    printf("%s-%d\n",insAllele[tti],insAlleleCnt[tti]);
  }

  //printf("%s",s[0]);
}
