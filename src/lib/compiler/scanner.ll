/* -*-C++-*- */
%{
/* Note the use of trailing context in the rules for the functions,
   allowing function names to be used as variable names.  This occurs
   in the "otree" example where "mean" is used as a variable name.
*/

#include <compiler/ParseTree.h>
#include "remap.h"
#include "parser.h"

%}

EXPONENT	[eE][+-][0-9]+
BRACKET		[ \t]*\(

%x COMMENT

%option nounput
%option yylineno

%%
"var"			return VAR;
"data"                  return DATA;
"model"			return MODEL;

"for"/{BRACKET}        	return FOR;
"in"/{BRACKET}		return IN;
"in"                    return IN;
";"			return ';';
","			return ',';
":"			return ':';
"["			return '[';
"]"			return ']';
"("			return '(';
")"			return ')';
"{"			return '{';
"}"			return '}';
"="			return '=';
"~"			return '~';
"<-"			return ARROW;
"+"			return '+';
"-"			return '-';
"*"			return '*';
"/"			return '/';

"T"/{BRACKET}           return 'T';

"/*"                    BEGIN(COMMENT);
<COMMENT>[^*]*          /* Eat up anything that's not a '*'   */
<COMMENT>"*"+[^*/n]*    /* Eat up '*'s not followed by a '/'  */
<COMMENT>"*"+"/"        BEGIN(INITIAL);

[ \t\r\n]+              /* Eat whitespace */
"#".*\n                 /* Eat comments */

([0-9]+){EXPONENT} {
  yylval.val = atof(yytext); return DOUBLE;
}
([0-9]+"."[0-9]*){EXPONENT} {
  yylval.val = atof(yytext); return DOUBLE;
}
([0-9]+"."[0-9]*) {
  yylval.val = atof(yytext); return DOUBLE;
}
("."[0-9]+){EXPONENT} {
  yylval.val = atof(yytext); return DOUBLE;
}
("."[0-9]+) {
  yylval.val = atof(yytext); return DOUBLE;
}
([0-9]+) {
  yylval.val = atof(yytext); return DOUBLE;
}


([a-zA-Z]+)/{BRACKET} {
  yylval.stringptr = new std::string(yytext);
  return FUNC;
}

[a-zA-Z]+[a-zA-Z0-9\.]*	 {
  yylval.stringptr = new std::string(yytext);
  return NAME;
}

<<EOF>>	yyterminate();

%%

int yywrap()
{
  return 1;
}
