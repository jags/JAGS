/* -*-C++-*- */
%{
#include <iostream>
#include <string>
#include <cstring>
#include <compiler/ParseTree.h>
#include "parser.hh"

    std::vector<FILE*> file_stack;
    int command_buffer_count = 0;
    int buffer_count = 0;
    void return_to_main_buffer();
    void close_buffer();
    std::string ExpandFileName(char const *name);
%}

%option prefix="zz"

EXPONENT	[eE][+-][0-9]+
BRACKET         [ \t]*\(

%s RDATA
%x COMMENT
%x SYSTEM

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
forceadapt			zzlval.intval=FORCEADAPT; return FORCEADAPT;
autoadapt			zzlval.intval=AUTOADAPT; return AUTOADAPT;
by                      zzlval.intval=BY; return BY;

monitor			zzlval.intval=MONITOR; return MONITOR;
monitors		zzlval.intval=MONITORS; return MONITORS;
nodenames		zzlval.intval=NODENAMES; return NODENAMES;
type			zzlval.intval=TYPE; return TYPE;
set			zzlval.intval=SET; return SET;
clear			zzlval.intval=CLEAR; return CLEAR;
thin			zzlval.intval=THIN; return THIN;
sampler                 zzlval.intval=SAMPLER; return SAMPLER;
samplers                zzlval.intval=SAMPLERS; return SAMPLERS;
rng                     zzlval.intval=RNGTOK; return RNGTOK;
factory                 zzlval.intval=FACTORY; return FACTORY;
factories               zzlval.intval=FACTORIES; return FACTORIES;
modules                 zzlval.intval=MODULES; return MODULES;
seed                    zzlval.intval=SEED; return SEED;

coda			zzlval.intval=CODA; return CODA;
stem			zzlval.intval=STEM; return STEM;

load                    zzlval.intval=LOAD; return LOAD;
unload                  zzlval.intval=UNLOAD; return UNLOAD;
exit			zzlval.intval=EXIT; return EXIT;

dir                     zzlval.intval=DIRECTORY; return DIRECTORY;
cd                      zzlval.intval=CD; return CD;
pwd                     zzlval.intval=PWD; return PWD;
run                     zzlval.intval=RUN; return RUN;

<INITIAL>"list"         zzlval.intval=LIST; return LIST;

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
<COMMENT>[^*]*          /* Eat up anything that's not a '*' */
<COMMENT>"*"+[^*/]*     /* Eat up '*'s not followed by a '/'  */
<COMMENT>"*"+"/"        BEGIN(INITIAL);

<RDATA,INITIAL>[ \t\r\f]+ /* Eat whitespace */
<RDATA,INITIAL>"#".*\n  /* Eat single-line comments */
<RDATA>[\n]+            /* Eat newlines */
<INITIAL>[\n]           return ENDCMD;

<INITIAL>"system"       BEGIN(SYSTEM);
<SYSTEM>.+ {
  zzlval.stringptr = new std::string(zztext);
  return SYSCMD;
}
<SYSTEM>[\n] {
    BEGIN(INITIAL); return ENDCMD;
 }

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
<RDATA>".Data"                  return DOTDATA;
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

_observations_ { 
  // Exception to the name rule for monitoring observed stochastic nodes:
  zzlval.stringptr = new std::string(zztext);
  return NAME;
}

[a-zA-Z0-9\._\-\\\/]+ {
  // More general than NAME which is limited to syntactically valid BUGS names
  zzlval.stringptr = new std::string(zztext);
  return STRING;
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
    if (command_buffer_count) {
	close_buffer();
	--command_buffer_count;
    }
    else {
	yyterminate();
    }
    return ENDSCRIPT;
}
<RDATA><<EOF>> {
    close_buffer();
    BEGIN(INITIAL);
    return ENDDATA;
}
%%

int zzwrap()
{
  return 1;
}

void push_file(FILE *file) {
    file_stack.push_back(file);
}

void pop_file() {
    if (file_stack.empty())
        return;

    fclose(file_stack.back());
    file_stack.pop_back();
}


bool open_data_buffer(std::string const *name) {
    FILE *file = fopen(ExpandFileName(name->c_str()).c_str(),"r");
    if (file) {
	zzpush_buffer_state(zz_create_buffer(file, YY_BUF_SIZE));
	push_file(file);
        ++buffer_count;
	BEGIN(RDATA);
	return true;
    }
    else {
	return false;
    }
}

void close_buffer() {
    zzpop_buffer_state();
    pop_file();
    --buffer_count;
}
    
bool open_command_buffer(std::string const *name) {
    FILE *file = fopen(name->c_str(),"r");
    if (file) {
	file_stack.push_back(file);
	zzpush_buffer_state(zz_create_buffer(file, YY_BUF_SIZE ));
	++command_buffer_count;
        ++buffer_count;
	return true;
    }
    else {
	return false;
    }
}

void return_to_main_buffer() {
    /* Clear all buffers in the stack and return to the first */
    while(buffer_count) {
	zzpop_buffer_state();
	--buffer_count;
    }
    command_buffer_count = 0;
    BEGIN(INITIAL);
}

