#include <pcre.h>
#include <stdio.h>


typedef struct {  // {{{
  pcre* reCompiled; 
  pcre_extra *pcreExtra;
} elmer_pcre_t;  // }}}


//TODO: emit messages according to msglevel
void pcre_f_compile(char* regStr, int msglevel, elmer_pcre_t* recompiled) { // {{{

  const char *pcreErrorStr;
  int pcreErrorOffset;

  recompiled->reCompiled = pcre_compile(regStr, 0, &pcreErrorStr, &pcreErrorOffset, NULL);
  recompiled->pcreExtra = pcre_study(recompiled->reCompiled, 0, &pcreErrorStr);
} // }}}


void pcre_f_match_compiled(char* Str, int StrLen, int* ovector, int maxmatch, elmer_pcre_t* compiled) { // {{{

  int pcreExecRet; 

  pcreExecRet = pcre_exec(compiled->reCompiled, NULL, Str, StrLen, 0, 0, ovector, 3*maxmatch);

} //  }}}
