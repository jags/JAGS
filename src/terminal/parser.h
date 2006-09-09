/* A Bison parser, made by GNU Bison 2.1.  */

/* Skeleton parser for Yacc-like parsing with Bison,
   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005 Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

/* As a special exception, when this file is copied by Bison into a
   Bison output file, you may use that output file without restriction.
   This special exception was added by the Free Software Foundation
   in version 1.24 of Bison.  */

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     INT = 258,
     DOUBLE = 259,
     NAME = 260,
     STRING = 261,
     ENDCMD = 262,
     MODEL = 263,
     DATA = 264,
     IN = 265,
     TO = 266,
     INITS = 267,
     PARAMETERS = 268,
     COMPILE = 269,
     INITIALIZE = 270,
     UPDATE = 271,
     BY = 272,
     MONITOR = 273,
     SET = 274,
     CLEAR = 275,
     THIN = 276,
     CODA = 277,
     STEM = 278,
     EXIT = 279,
     NCHAINS = 280,
     CHAIN = 281,
     LOAD = 282,
     LIST = 283,
     STRUCTURE = 284,
     DIM = 285,
     NA = 286,
     R_NULL = 287,
     DIMNAMES = 288,
     ITER = 289,
     ARROW = 290,
     ENDDATA = 291,
     ASINTEGER = 292
   };
#endif
/* Tokens.  */
#define INT 258
#define DOUBLE 259
#define NAME 260
#define STRING 261
#define ENDCMD 262
#define MODEL 263
#define DATA 264
#define IN 265
#define TO 266
#define INITS 267
#define PARAMETERS 268
#define COMPILE 269
#define INITIALIZE 270
#define UPDATE 271
#define BY 272
#define MONITOR 273
#define SET 274
#define CLEAR 275
#define THIN 276
#define CODA 277
#define STEM 278
#define EXIT 279
#define NCHAINS 280
#define CHAIN 281
#define LOAD 282
#define LIST 283
#define STRUCTURE 284
#define DIM 285
#define NA 286
#define R_NULL 287
#define DIMNAMES 288
#define ITER 289
#define ARROW 290
#define ENDDATA 291
#define ASINTEGER 292




#if ! defined (YYSTYPE) && ! defined (YYSTYPE_IS_DECLARED)
#line 54 "parser.yy"
typedef union YYSTYPE {
  int intval;
  double val;
  std::string *stringptr;
  ParseTree *ptree;
  std::vector<ParseTree*> *pvec;
  std::vector<double> *vec;
  std::vector<long> *ivec;
} YYSTYPE;
/* Line 1447 of yacc.c.  */
#line 122 "parser.h"
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE zzlval;



