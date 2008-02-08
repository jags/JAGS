/* -*-C++-*- */
%{
#include <iostream>
#include <cstring>
#include <compiler/ParseTree.h>
#include "parser.h"
    
    std::vector<YY_BUFFER_STATE> buffer_stack;
    std::vector<bool> int_stack;
    void push_buffer(YY_BUFFER_STATE);
    void pop_buffer();
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
adapt			zzlval.intval=ADAPT; return ADAPT;
by                      zzlval.intval=BY; return BY;

monitor			zzlval.intval=MONITOR; return MONITOR;
monitors		zzlval.intval=MONITORS; return MONITORS;
type			zzlval.intval=TYPE; return TYPE;
set			zzlval.intval=SET; return SET;
clear			zzlval.intval=CLEAR; return CLEAR;
thin			zzlval.intval=THIN; return THIN;
samplers                zzlval.intval=SAMPLERS; return SAMPLERS;

coda			zzlval.intval=CODA; return CODA;
stem			zzlval.intval=STEM; return STEM;

load                    zzlval.intval=LOAD; return LOAD;
exit			zzlval.intval=EXIT; return EXIT;

dir                     zzlval.intval=DIRECTORY; return DIRECTORY;
cd                      zzlval.intval=CD; return CD;
pwd                     zzlval.intval=PWD; return PWD;
run                     zzlval.intval=RUN; return RUN;

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
"`"                     return '`';

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
<RDATA>"-"?[0-9]+"L" {
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

\'[^\']*\' {
    zzlval.stringptr = new std::string(zztext);
    zzlval.stringptr->erase(zzlval.stringptr->size() - 1, 1);
    zzlval.stringptr->erase(0,1);
    return STRING;
}

<INITIAL><<EOF>> {
    if (buffer_stack.empty()) {
	yyterminate();
    }
    else {
	pop_buffer();
    }
    return ENDSCRIPT;
}
<RDATA><<EOF>> {
    pop_buffer();
    BEGIN(INITIAL);
    return ENDDATA;
}
%%

int zzwrap()
{
  return 1;
}

bool open_data_buffer(std::string const *name) {
    FILE *file = fopen(name->c_str(),"r");
    if (file) {
	push_buffer(YY_CURRENT_BUFFER);
	yy_switch_to_buffer(yy_create_buffer(file, YY_BUF_SIZE ) );
	BEGIN(RDATA);
	return true;
    }
    else {
	return false;
    }
}

bool open_command_buffer(std::string const *name) {
    FILE *file = fopen(name->c_str(),"r");
    if (file) {
	push_buffer(YY_CURRENT_BUFFER);
	yy_switch_to_buffer(yy_create_buffer(file, YY_BUF_SIZE ));
	return true;
    }
    else {
	return false;
    }
}

void push_buffer(YY_BUFFER_STATE new_buffer) {
    buffer_stack.push_back(new_buffer);
}

void pop_buffer() {
    if (!buffer_stack.empty()) {
	yy_delete_buffer (YY_CURRENT_BUFFER);
	yy_switch_to_buffer(buffer_stack.back());
	buffer_stack.pop_back();
    }
}

void return_to_main_buffer() {
    /* Clear all buffers in the stack and return to the first */
    if (!buffer_stack.empty()) {
	for (unsigned int i = 1; i < buffer_stack.size(); ++i) {
	    yy_delete_buffer (buffer_stack[i]);
	}
	yy_switch_to_buffer(buffer_stack[0]);
	buffer_stack.clear();
    }
    BEGIN(INITIAL);
}

