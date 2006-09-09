/* -*-C++-*- */
%{
#include <iostream>
#include <cstring>
#include <compiler/ParseTree.h>
#include "parser.h"
  YY_BUFFER_STATE main_buffer = 0;
  void return_to_main_buffer();
%}

%option prefix="zz"

EXPONENT	[eE][+-][0-9]+
BRACKET         [ \t]*\(

%s RDATA
%x COMMENT

%%

model                   zzlval.intval=MODEL; return MODEL;
data                    zzlval.intval=DATA; return DATA;
in                      zzlval.intval=IN; return IN;
to			zzlval.intval=TO; return TO;
inits                   zzlval.intval=INITS; return INITS;
parameters		zzlval.intval=PARAMETERS; return PARAMETERS;

compile			zzlval.intval=COMPILE; return COMPILE;
nchains                 zzlval.intval=NCHAINS; return NCHAINS;
initialize              zzlval.intval=INITIALIZE; return INITIALIZE;
chain                   zzlval.intval=CHAIN; return CHAIN;

update			zzlval.intval=UPDATE; return UPDATE;
by                      zzlval.intval=BY; return BY;

monitor			zzlval.intval=MONITOR; return MONITOR;
set			zzlval.intval=SET; return SET;
clear			zzlval.intval=CLEAR; return CLEAR;
thin			zzlval.intval=THIN; return THIN;

coda			zzlval.intval=CODA; return CODA;
stem			zzlval.intval=STEM; return STEM;

load                    zzlval.intval=LOAD; return LOAD;
exit			zzlval.intval=EXIT; return EXIT;

","			return ',';
"["			return '[';
"]"			return ']';
"("			return '(';
")"			return ')';
"="                     return '=';
"<-"			return ARROW;
"*"			return '*';
":"                     return ':';
";"                     return ';';

"/*"                    BEGIN(COMMENT);
<COMMENT>[^*]*          /* Eat up anything that's not a '*'   */
<COMMENT>"*"+[^*/n]*    /* Eat up '*'s not followed by a '/'  */
<COMMENT>"*"+"/"        BEGIN(INITIAL);

<RDATA,INITIAL>[ \t\r]+ /* Eat whitespace */
<RDATA,INITIAL>"#".*\n          /* Eat comments */
<RDATA>[\n]+            /* Eat newlines */
<INITIAL>[\n]           return ENDCMD;

<INITIAL,RDATA>"-"?([0-9]+){EXPONENT}  {
  zzlval.val = atof(zztext); return DOUBLE;
}
<INITIAL,RDATA>"-"?([0-9]+"."[0-9]*){EXPONENT}  {
  zzlval.val = atof(zztext); return DOUBLE;
}
<INITIAL,RDATA>"-"?([0-9]+"."[0-9]*)  {
  zzlval.val = atof(zztext); return DOUBLE;
}
<INITIAL,RDATA>"-"?("."[0-9]+){EXPONENT}  {
  zzlval.val = atof(zztext); return DOUBLE;
}
<INITIAL,RDATA>"-"?("."[0-9]+)  {
  zzlval.val = atof(zztext); return DOUBLE;
}
<RDATA>"-"?[0-9]+  {
  zzlval.val = atof(zztext); return DOUBLE;
}
<INITIAL>"-"?[0-9]+	{
  zzlval.intval = atoi(zztext); return INT;
}

<RDATA>".Dim"                   return DIM;
<RDATA>"NA"		        return NA;
<RDATA>"NULL"		        return R_NULL;
<RDATA>"structure"/{BRACKET}	return STRUCTURE;
<RDATA>"list"/{BRACKET}		return LIST;
<RDATA>"c"/{BRACKET} 		return 'c';
<RDATA>"as.integer"/{BRACKET}   return ASINTEGER;

[a-zA-Z\.]+[a-zA-Z0-9\._\-\\\/]* { 
  zzlval.stringptr = new std::string(zztext);
  return NAME;
}

\"[^\"]*\"  {
  zzlval.stringptr = new std::string(zztext);
  // Remove enclosing quotes
  zzlval.stringptr->erase(zzlval.stringptr->size() - 1, 1);
  zzlval.stringptr->erase(0, 1);
  return STRING;
}

<INITIAL><<EOF>>		yyterminate();
<RDATA><<EOF>> {
  return_to_main_buffer();
  return ENDDATA;
}
%%

int zzwrap()
{
  return 1;
}

bool open_data_buffer(std::string const *name) {
  if (main_buffer) {
    return false;
  }
  FILE *file = fopen(name->c_str(),"r");
  if (file) {
    main_buffer = YY_CURRENT_BUFFER;
    yy_switch_to_buffer(yy_create_buffer(file, YY_BUF_SIZE ) );
    BEGIN(RDATA);
    return true;
  }
  else {
    return false;
  }
}

void return_to_main_buffer() {
  if (main_buffer && (main_buffer != YY_CURRENT_BUFFER)) {
    yy_delete_buffer (YY_CURRENT_BUFFER);
    yy_switch_to_buffer (main_buffer);
    main_buffer = 0;
  }
  BEGIN(INITIAL);  
}
