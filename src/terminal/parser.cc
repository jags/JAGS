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

/* Written by Richard Stallman by simplifying the original so called
   ``semantic'' parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "2.1"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Using locations.  */
#define YYLSP_NEEDED 0

/* Substitute the variable and function names.  */
#define yyparse zzparse
#define yylex   zzlex
#define yyerror zzerror
#define yylval  zzlval
#define yychar  zzchar
#define yydebug zzdebug
#define yynerrs zznerrs


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




/* Copy the first part of user declarations.  */
#line 2 "parser.yy"

#include <config.h>

/* #include <Rmath.h> */

#include <cstdio>
#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <fstream>

#include <Console.h>
#include <compiler/ParseTree.h>
#include <sarray/SArray.h>
#include <cstdlib>
#include <cstring>
#include <ltdl.h>

#include "ReadData.h"
//#include "JRNG.h"

    typedef void(*pt2Func)();

    int zzerror(const char *);
    int zzlex();
#define YYERROR_VERBOSE 0
    static Console *console;
    static bool interactive = true;
    void setName(ParseTree *p, std::string *name);
    std::map<std::string, SArray> _data_table;
    bool open_data_buffer(std::string const *name);
    void return_to_main_buffer();
    void setMonitor(ParseTree const *var, int thin);
    void clearMonitor(ParseTree const *var);
    void doCoda (ParseTree const *var, std::string const &stem);
    void doAllCoda (std::string const &stem);
    void doDump (std::string const &file, DumpType type, unsigned int chain);

    static void errordump();
    static void updatestar(long niter, long refresh, int width);
    static void setParameters(ParseTree *p, ParseTree *param1);
    static void setParameters(ParseTree *p, std::vector<ParseTree*> *parameters);
    static void setParameters(ParseTree *p, ParseTree *param1, ParseTree *param2);
    static void loadModule(std::string const &name);
    

/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* Enabling the token table.  */
#ifndef YYTOKEN_TABLE
# define YYTOKEN_TABLE 0
#endif

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
/* Line 196 of yacc.c.  */
#line 226 "parser.cc"
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 219 of yacc.c.  */
#line 238 "parser.cc"

#if ! defined (YYSIZE_T) && defined (__SIZE_TYPE__)
# define YYSIZE_T __SIZE_TYPE__
#endif
#if ! defined (YYSIZE_T) && defined (size_t)
# define YYSIZE_T size_t
#endif
#if ! defined (YYSIZE_T) && (defined (__STDC__) || defined (__cplusplus))
# include <stddef.h> /* INFRINGES ON USER NAME SPACE */
# define YYSIZE_T size_t
#endif
#if ! defined (YYSIZE_T)
# define YYSIZE_T unsigned int
#endif

#ifndef YY_
# if YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

#if ! defined (yyoverflow) || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if defined (__STDC__) || defined (__cplusplus)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#     define YYINCLUDED_STDLIB_H
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning. */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2005 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM ((YYSIZE_T) -1)
#  endif
#  ifdef __cplusplus
extern "C" {
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if (! defined (malloc) && ! defined (YYINCLUDED_STDLIB_H) \
	&& (defined (__STDC__) || defined (__cplusplus)))
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if (! defined (free) && ! defined (YYINCLUDED_STDLIB_H) \
	&& (defined (__STDC__) || defined (__cplusplus)))
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifdef __cplusplus
}
#  endif
# endif
#endif /* ! defined (yyoverflow) || YYERROR_VERBOSE */


#if (! defined (yyoverflow) \
     && (! defined (__cplusplus) \
	 || (defined (YYSTYPE_IS_TRIVIAL) && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  short int yyss;
  YYSTYPE yyvs;
  };

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (short int) + sizeof (YYSTYPE))			\
      + YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined (__GNUC__) && 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  YYSIZE_T yyi;				\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (0)
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack)					\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack, Stack, yysize);				\
	Stack = &yyptr->Stack;						\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (0)

#endif

#if defined (__STDC__) || defined (__cplusplus)
   typedef signed char yysigned_char;
#else
   typedef short int yysigned_char;
#endif

/* YYFINAL -- State number of the termination state. */
#define YYFINAL  2
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   180

/* YYNTOKENS -- Number of terminals. */
#define YYNTOKENS  48
/* YYNNTS -- Number of nonterminals. */
#define YYNNTS  44
/* YYNRULES -- Number of rules. */
#define YYNRULES  99
/* YYNRULES -- Number of states. */
#define YYNSTATES  196

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   292

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const unsigned char yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      39,    40,    44,     2,    38,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    43,    45,
       2,    46,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    41,     2,    42,     2,     2,     2,     2,     2,    47,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const unsigned short int yyprhs[] =
{
       0,     0,     3,     4,     7,     9,    12,    15,    17,    19,
      21,    23,    25,    27,    29,    31,    33,    35,    37,    39,
      41,    45,    49,    51,    55,    59,    62,    66,    75,    77,
      81,    90,    94,    98,   100,   107,   109,   112,   120,   122,
     124,   129,   131,   135,   137,   141,   143,   145,   147,   151,
     160,   163,   171,   175,   177,   179,   182,   190,   193,   201,
     204,   206,   209,   213,   217,   221,   225,   232,   234,   236,
     240,   244,   248,   250,   252,   257,   259,   264,   266,   270,
     272,   274,   278,   280,   284,   286,   290,   292,   297,   299,
     304,   309,   311,   313,   318,   320,   324,   326,   331,   333
};

/* YYRHS -- A `-1'-separated list of the rules' RHS. */
static const yysigned_char yyrhs[] =
{
      49,     0,    -1,    -1,    49,    50,    -1,     7,    -1,    51,
       7,    -1,     1,     7,    -1,    52,    -1,    53,    -1,    54,
      -1,    56,    -1,    57,    -1,    58,    -1,    60,    -1,    61,
      -1,    62,    -1,    68,    -1,    72,    -1,    73,    -1,    63,
      -1,     8,    10,    71,    -1,    55,    74,    36,    -1,    55,
      -1,     9,    11,    71,    -1,     9,    10,    71,    -1,     9,
      20,    -1,    59,    74,    36,    -1,    59,    74,    36,    38,
      26,    39,     3,    40,    -1,    59,    -1,    13,    11,    71,
      -1,    13,    11,    71,    38,    26,    39,     3,    40,    -1,
      13,    10,    71,    -1,    12,    10,    71,    -1,    14,    -1,
      14,    38,    25,    39,     3,    40,    -1,    15,    -1,    16,
       3,    -1,    16,     3,    38,    17,    39,     3,    40,    -1,
      24,    -1,     5,    -1,     5,    41,    65,    42,    -1,    66,
      -1,    65,    38,    66,    -1,    67,    -1,    67,    43,    67,
      -1,     3,    -1,    69,    -1,    70,    -1,    18,    19,    64,
      -1,    18,    19,    64,    38,    21,    39,     3,    40,    -1,
      18,    64,    -1,    18,    64,    38,    21,    39,     3,    40,
      -1,    18,    20,    64,    -1,     5,    -1,     6,    -1,    22,
      64,    -1,    22,    64,    38,    23,    39,    71,    40,    -1,
      22,    44,    -1,    22,    44,    38,    23,    39,    71,    40,
      -1,    27,    71,    -1,    75,    -1,    74,    75,    -1,    74,
      45,    75,    -1,     6,    35,    76,    -1,     6,    35,    79,
      -1,     6,    35,     6,    -1,    29,    39,    79,    38,    77,
      40,    -1,    78,    -1,    84,    -1,    77,    38,    84,    -1,
      77,    38,    78,    -1,    30,    46,    79,    -1,    80,    -1,
      81,    -1,    37,    39,    81,    40,    -1,    83,    -1,    47,
      39,    82,    40,    -1,    83,    -1,    82,    38,    83,    -1,
       4,    -1,    31,    -1,     5,    46,    87,    -1,    86,    -1,
      85,    38,    86,    -1,    87,    -1,     5,    46,    87,    -1,
      88,    -1,    37,    39,    88,    40,    -1,    90,    -1,    28,
      39,    85,    40,    -1,    29,    39,    85,    40,    -1,    32,
      -1,     4,    -1,    47,    39,    89,    40,    -1,     4,    -1,
      89,    38,     4,    -1,     6,    -1,    47,    39,    91,    40,
      -1,     6,    -1,    91,    38,     6,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const unsigned short int yyrline[] =
{
       0,   114,   114,   115,   118,   119,   120,   123,   124,   125,
     126,   127,   128,   129,   130,   131,   132,   133,   134,   135,
     138,   151,   159,   162,   168,   179,   185,   198,   211,   214,
     218,   224,   233,   247,   250,   255,   262,   267,   272,   275,
     278,   284,   287,   292,   295,   301,   304,   305,   308,   311,
     314,   317,   322,   331,   332,   335,   338,   341,   344,   349,
     354,   357,   360,   365,   368,   373,   381,   392,   393,   394,
     395,   398,   403,   404,   407,   410,   414,   420,   421,   424,
     425,   431,   434,   435,   438,   439,   442,   443,   444,   445,
     446,   447,   450,   451,   454,   455,   458,   459,   462,   463
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals. */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "INT", "DOUBLE", "NAME", "STRING",
  "ENDCMD", "MODEL", "DATA", "IN", "TO", "INITS", "PARAMETERS", "COMPILE",
  "INITIALIZE", "UPDATE", "BY", "MONITOR", "SET", "CLEAR", "THIN", "CODA",
  "STEM", "EXIT", "NCHAINS", "CHAIN", "LOAD", "LIST", "STRUCTURE", "DIM",
  "NA", "R_NULL", "DIMNAMES", "ITER", "ARROW", "ENDDATA", "ASINTEGER",
  "','", "'('", "')'", "'['", "']'", "':'", "'*'", "';'", "'='", "'c'",
  "$accept", "input", "line", "command", "model", "data_in", "data_to",
  "data", "data_clear", "parameters_in", "parameters_to", "parameters",
  "compile", "initialize", "update", "exit", "var", "range_list",
  "range_element", "index", "monitor", "monitor_set", "monitor_clear",
  "file_name", "coda", "load", "r_assignment_list", "r_assignment",
  "r_structure", "r_attribute_list", "r_dim", "r_collection",
  "r_integer_collection", "r_value_collection", "r_value_list", "r_value",
  "r_generic_attribute", "r_generic_list", "r_generic_list_element",
  "r_generic_vector", "r_numeric_vector", "r_double_list",
  "r_character_vector", "r_string_list", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const unsigned short int yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,    44,    40,
      41,    91,    93,    58,    42,    59,    61,    99
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const unsigned char yyr1[] =
{
       0,    48,    49,    49,    50,    50,    50,    51,    51,    51,
      51,    51,    51,    51,    51,    51,    51,    51,    51,    51,
      52,    53,    53,    54,    55,    56,    57,    57,    57,    58,
      58,    59,    59,    60,    60,    61,    62,    62,    63,    64,
      64,    65,    65,    66,    66,    67,    68,    68,    69,    69,
      69,    69,    70,    71,    71,    72,    72,    72,    72,    73,
      74,    74,    74,    75,    75,    75,    76,    77,    77,    77,
      77,    78,    79,    79,    80,    81,    81,    82,    82,    83,
      83,    84,    85,    85,    86,    86,    87,    87,    87,    87,
      87,    87,    88,    88,    89,    89,    90,    90,    91,    91
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const unsigned char yyr2[] =
{
       0,     2,     0,     2,     1,     2,     2,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       3,     3,     1,     3,     3,     2,     3,     8,     1,     3,
       8,     3,     3,     1,     6,     1,     2,     7,     1,     1,
       4,     1,     3,     1,     3,     1,     1,     1,     3,     8,
       2,     7,     3,     1,     1,     2,     7,     2,     7,     2,
       1,     2,     3,     3,     3,     3,     6,     1,     1,     3,
       3,     3,     1,     1,     4,     1,     4,     1,     3,     1,
       1,     3,     1,     3,     1,     3,     1,     4,     1,     4,
       4,     1,     1,     4,     1,     3,     1,     4,     1,     3
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const unsigned char yydefact[] =
{
       2,     0,     1,     0,     4,     0,     0,     0,     0,    33,
      35,     0,     0,     0,    38,     0,     3,     0,     7,     8,
       9,    22,    10,    11,    12,    28,    13,    14,    15,    19,
      16,    46,    47,    17,    18,     6,     0,     0,     0,    25,
       0,     0,     0,     0,    36,    39,     0,     0,    50,    57,
      55,    53,    54,    59,     5,     0,     0,    60,     0,    20,
      24,    23,    32,    31,    29,     0,     0,     0,    48,    52,
       0,     0,     0,     0,    21,     0,    61,    26,     0,     0,
       0,    45,     0,    41,    43,     0,     0,     0,     0,    79,
      65,     0,    80,     0,     0,    63,    64,    72,    73,    75,
      62,     0,     0,     0,     0,     0,    40,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    34,     0,    42,
      44,     0,     0,     0,     0,     0,     0,     0,    77,     0,
       0,    37,     0,    51,    58,    56,     0,    74,     0,    76,
       0,    30,    49,     0,     0,     0,    67,    68,    78,    27,
       0,     0,     0,    66,    92,    96,     0,     0,    91,     0,
       0,    81,    86,    88,    71,    70,    69,     0,     0,     0,
       0,     0,     0,    82,    84,     0,     0,     0,    94,    98,
       0,     0,     0,     0,    89,    90,     0,    87,     0,    93,
       0,    97,    85,    83,    95,    99
};

/* YYDEFGOTO[NTERM-NUM]. */
static const short int yydefgoto[] =
{
      -1,     1,    16,    17,    18,    19,    20,    21,    22,    23,
      24,    25,    26,    27,    28,    29,    48,    82,    83,    84,
      30,    31,    32,    53,    33,    34,    56,    57,    95,   145,
     146,    96,    97,    98,   127,    99,   147,   172,   173,   174,
     162,   180,   163,   181
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -126
static const short int yypact[] =
{
    -126,    69,  -126,     4,  -126,     7,    78,    19,    97,    28,
    -126,    56,    60,    18,  -126,    67,  -126,    79,  -126,  -126,
    -126,    88,  -126,  -126,  -126,    88,  -126,  -126,  -126,  -126,
    -126,  -126,  -126,  -126,  -126,  -126,    67,    67,    67,  -126,
      67,    67,    67,    84,    72,    70,   107,   107,    75,    76,
      77,  -126,  -126,  -126,  -126,    81,    -3,  -126,    13,  -126,
    -126,  -126,  -126,  -126,    80,    82,   100,   116,    85,  -126,
      99,   101,   102,    14,  -126,    88,  -126,    89,    96,   123,
      90,  -126,   -16,  -126,    87,   110,    93,    94,    95,  -126,
    -126,    98,  -126,   103,   104,  -126,  -126,  -126,  -126,  -126,
    -126,   109,   105,   106,   125,   116,  -126,   116,   108,   133,
      67,    67,    17,     8,    40,   111,   135,  -126,   112,  -126,
    -126,   136,   113,   114,   115,   118,   117,   -24,  -126,   137,
     119,  -126,   120,  -126,  -126,  -126,    22,  -126,    40,  -126,
     121,  -126,  -126,   122,   124,    52,  -126,  -126,  -126,  -126,
       9,    17,    22,  -126,  -126,  -126,   126,   127,  -126,   128,
     130,  -126,  -126,  -126,  -126,  -126,  -126,     3,     3,     6,
      91,   129,    61,  -126,  -126,    62,   132,   134,  -126,  -126,
      65,    66,     9,     3,  -126,  -126,   141,  -126,   144,  -126,
     143,  -126,  -126,  -126,  -126,  -126
};

/* YYPGOTO[NTERM-NUM].  */
static const short int yypgoto[] =
{
    -126,  -126,  -126,  -126,  -126,  -126,  -126,  -126,  -126,  -126,
    -126,  -126,  -126,  -126,  -126,  -126,    21,  -126,    36,    44,
    -126,  -126,  -126,   -36,  -126,  -126,   138,   -28,  -126,  -126,
      10,   -88,  -126,    45,  -126,   -78,    12,     5,   -11,  -125,
      11,  -126,  -126,  -126
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -1
static const unsigned char yytable[] =
{
      59,    60,    61,    55,    62,    63,    64,   154,   171,   155,
     154,    35,    89,   154,   138,   155,   139,    36,    89,    55,
      90,    89,   105,    45,   125,   161,   106,   143,    76,    40,
      76,   156,   157,    74,    50,   158,   128,   156,   157,    92,
     159,   158,    75,    91,    89,    92,   159,   100,    92,    77,
     160,    93,   144,   176,    93,    94,   160,   192,    75,    44,
     148,    94,    49,   164,    94,    45,    43,    68,    69,     2,
       3,    92,    51,    52,   123,   124,     4,     5,     6,    46,
      47,     7,     8,     9,    10,    11,    54,    12,    37,    38,
     152,    13,   153,    14,    55,   178,    15,   179,    39,   183,
     183,   184,   185,   188,   190,   189,   191,    41,    42,    65,
      66,    67,    45,    70,    71,    72,    73,    80,    78,    81,
      86,    79,   102,    85,    87,    88,   103,   101,   118,   104,
     107,   108,   109,   110,   111,   115,   122,   112,   130,   132,
     140,   119,   113,   114,   116,   178,   117,   121,   194,   195,
     129,   120,   131,   133,   134,   135,   136,   137,   126,   141,
     142,   149,   165,    58,   166,   167,   168,   169,   150,   170,
     151,   186,   193,   175,   187,   182,     0,     0,     0,     0,
     177
};

static const short int yycheck[] =
{
      36,    37,    38,     6,    40,    41,    42,     4,     5,     6,
       4,     7,     4,     4,    38,     6,    40,    10,     4,     6,
       6,     4,    38,     5,   112,   150,    42,     5,    56,    10,
      58,    28,    29,    36,    13,    32,   114,    28,    29,    31,
      37,    32,    45,    29,     4,    31,    37,    75,    31,    36,
      47,    37,    30,    47,    37,    47,    47,   182,    45,     3,
     138,    47,    44,   151,    47,     5,    38,    46,    47,     0,
       1,    31,     5,     6,   110,   111,     7,     8,     9,    19,
      20,    12,    13,    14,    15,    16,     7,    18,    10,    11,
      38,    22,    40,    24,     6,     4,    27,     6,    20,    38,
      38,    40,    40,    38,    38,    40,    40,    10,    11,    25,
      38,    41,     5,    38,    38,    38,    35,    17,    38,     3,
      21,    39,    26,    38,    23,    23,     3,    38,     3,    39,
      43,    21,    39,    39,    39,    26,     3,    39,     3,     3,
       3,   105,    39,    39,    39,     4,    40,    39,     4,     6,
      39,   107,    40,    40,    40,    40,    38,    40,   113,    40,
      40,    40,   152,    25,   152,    39,    39,    39,    46,    39,
      46,    39,   183,   168,    40,    46,    -1,    -1,    -1,    -1,
     169
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const unsigned char yystos[] =
{
       0,    49,     0,     1,     7,     8,     9,    12,    13,    14,
      15,    16,    18,    22,    24,    27,    50,    51,    52,    53,
      54,    55,    56,    57,    58,    59,    60,    61,    62,    63,
      68,    69,    70,    72,    73,     7,    10,    10,    11,    20,
      10,    10,    11,    38,     3,     5,    19,    20,    64,    44,
      64,     5,     6,    71,     7,     6,    74,    75,    74,    71,
      71,    71,    71,    71,    71,    25,    38,    41,    64,    64,
      38,    38,    38,    35,    36,    45,    75,    36,    38,    39,
      17,     3,    65,    66,    67,    38,    21,    23,    23,     4,
       6,    29,    31,    37,    47,    76,    79,    80,    81,    83,
      75,    38,    26,     3,    39,    38,    42,    43,    21,    39,
      39,    39,    39,    39,    39,    26,    39,    40,     3,    66,
      67,    39,     3,    71,    71,    79,    81,    82,    83,    39,
       3,    40,     3,    40,    40,    40,    38,    40,    38,    40,
       3,    40,    40,     5,    30,    77,    78,    84,    83,    40,
      46,    46,    38,    40,     4,     6,    28,    29,    32,    37,
      47,    87,    88,    90,    79,    78,    84,    39,    39,    39,
      39,     5,    85,    86,    87,    85,    47,    88,     4,     6,
      89,    91,    46,    38,    40,    40,    39,    40,    38,    40,
      38,    40,    87,    86,     4,     6
};

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */

#define YYFAIL		goto yyerrlab

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)					\
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    {								\
      yychar = (Token);						\
      yylval = (Value);						\
      yytoken = YYTRANSLATE (yychar);				\
      YYPOPSTACK;						\
      goto yybackup;						\
    }								\
  else								\
    {								\
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (0)


#define YYTERROR	1
#define YYERRCODE	256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#define YYRHSLOC(Rhs, K) ((Rhs)[K])
#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)				\
    do									\
      if (N)								\
	{								\
	  (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;	\
	  (Current).first_column = YYRHSLOC (Rhs, 1).first_column;	\
	  (Current).last_line    = YYRHSLOC (Rhs, N).last_line;		\
	  (Current).last_column  = YYRHSLOC (Rhs, N).last_column;	\
	}								\
      else								\
	{								\
	  (Current).first_line   = (Current).last_line   =		\
	    YYRHSLOC (Rhs, 0).last_line;				\
	  (Current).first_column = (Current).last_column =		\
	    YYRHSLOC (Rhs, 0).last_column;				\
	}								\
    while (0)
#endif


/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
# if YYLTYPE_IS_TRIVIAL
#  define YY_LOCATION_PRINT(File, Loc)			\
     fprintf (File, "%d.%d-%d.%d",			\
              (Loc).first_line, (Loc).first_column,	\
              (Loc).last_line,  (Loc).last_column)
# else
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX yylex (YYLEX_PARAM)
#else
# define YYLEX yylex ()
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (0)

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)		\
do {								\
  if (yydebug)							\
    {								\
      YYFPRINTF (stderr, "%s ", Title);				\
      yysymprint (stderr,					\
                  Type, Value);	\
      YYFPRINTF (stderr, "\n");					\
    }								\
} while (0)

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yy_stack_print (short int *bottom, short int *top)
#else
static void
yy_stack_print (bottom, top)
    short int *bottom;
    short int *top;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (/* Nothing. */; bottom <= top; ++bottom)
    YYFPRINTF (stderr, " %d", *bottom);
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yy_reduce_print (int yyrule)
#else
static void
yy_reduce_print (yyrule)
    int yyrule;
#endif
{
  int yyi;
  unsigned long int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu), ",
             yyrule - 1, yylno);
  /* Print the symbols being reduced, and their result.  */
  for (yyi = yyprhs[yyrule]; 0 <= yyrhs[yyi]; yyi++)
    YYFPRINTF (stderr, "%s ", yytname[yyrhs[yyi]]);
  YYFPRINTF (stderr, "-> %s\n", yytname[yyr1[yyrule]]);
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (Rule);		\
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif



#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined (__GLIBC__) && defined (_STRING_H)
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
static YYSIZE_T
#   if defined (__STDC__) || defined (__cplusplus)
yystrlen (const char *yystr)
#   else
yystrlen (yystr)
     const char *yystr;
#   endif
{
  const char *yys = yystr;

  while (*yys++ != '\0')
    continue;

  return yys - yystr - 1;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined (__GLIBC__) && defined (_STRING_H) && defined (_GNU_SOURCE)
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
static char *
#   if defined (__STDC__) || defined (__cplusplus)
yystpcpy (char *yydest, const char *yysrc)
#   else
yystpcpy (yydest, yysrc)
     char *yydest;
     const char *yysrc;
#   endif
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      size_t yyn = 0;
      char const *yyp = yystr;

      for (;;)
	switch (*++yyp)
	  {
	  case '\'':
	  case ',':
	    goto do_not_strip_quotes;

	  case '\\':
	    if (*++yyp != '\\')
	      goto do_not_strip_quotes;
	    /* Fall through.  */
	  default:
	    if (yyres)
	      yyres[yyn] = *yyp;
	    yyn++;
	    break;

	  case '"':
	    if (yyres)
	      yyres[yyn] = '\0';
	    return yyn;
	  }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

#endif /* YYERROR_VERBOSE */



#if YYDEBUG
/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yysymprint (FILE *yyoutput, int yytype, YYSTYPE *yyvaluep)
#else
static void
yysymprint (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  /* Pacify ``unused variable'' warnings.  */
  (void) yyvaluep;

  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);


# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# endif
  switch (yytype)
    {
      default:
        break;
    }
  YYFPRINTF (yyoutput, ")");
}

#endif /* ! YYDEBUG */
/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yymsg, yytype, yyvaluep)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  /* Pacify ``unused variable'' warnings.  */
  (void) yyvaluep;

  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  switch (yytype)
    {

      default:
        break;
    }
}


/* Prevent warnings from -Wmissing-prototypes.  */

#ifdef YYPARSE_PARAM
# if defined (__STDC__) || defined (__cplusplus)
int yyparse (void *YYPARSE_PARAM);
# else
int yyparse ();
# endif
#else /* ! YYPARSE_PARAM */
#if defined (__STDC__) || defined (__cplusplus)
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */



/* The look-ahead symbol.  */
int yychar;

/* The semantic value of the look-ahead symbol.  */
YYSTYPE yylval;

/* Number of syntax errors so far.  */
int yynerrs;



/*----------.
| yyparse.  |
`----------*/

#ifdef YYPARSE_PARAM
# if defined (__STDC__) || defined (__cplusplus)
int yyparse (void *YYPARSE_PARAM)
# else
int yyparse (YYPARSE_PARAM)
  void *YYPARSE_PARAM;
# endif
#else /* ! YYPARSE_PARAM */
#if defined (__STDC__) || defined (__cplusplus)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{
  
  int yystate;
  int yyn;
  int yyresult;
  /* Number of tokens to shift before error messages enabled.  */
  int yyerrstatus;
  /* Look-ahead token as an internal (translated) token number.  */
  int yytoken = 0;

  /* Three stacks and their tools:
     `yyss': related to states,
     `yyvs': related to semantic values,
     `yyls': related to locations.

     Refer to the stacks thru separate pointers, to allow yyoverflow
     to reallocate them elsewhere.  */

  /* The state stack.  */
  short int yyssa[YYINITDEPTH];
  short int *yyss = yyssa;
  short int *yyssp;

  /* The semantic value stack.  */
  YYSTYPE yyvsa[YYINITDEPTH];
  YYSTYPE *yyvs = yyvsa;
  YYSTYPE *yyvsp;



#define YYPOPSTACK   (yyvsp--, yyssp--)

  YYSIZE_T yystacksize = YYINITDEPTH;

  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;


  /* When reducing, the number of symbols on the RHS of the reduced
     rule.  */
  int yylen;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY;		/* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  yyssp = yyss;
  yyvsp = yyvs;

  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed. so pushing a state here evens the stacks.
     */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack. Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	short int *yyss1 = yyss;


	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow (YY_("memory exhausted"),
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),

		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	short int *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyexhaustedlab;
	YYSTACK_RELOCATE (yyss);
	YYSTACK_RELOCATE (yyvs);

#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;


      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

/* Do appropriate processing given the current state.  */
/* Read a look-ahead token if we need one and don't already have one.  */
/* yyresume: */

  /* First try to decide what to do without reference to look-ahead token.  */

  yyn = yypact[yystate];
  if (yyn == YYPACT_NINF)
    goto yydefault;

  /* Not known => get a look-ahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid look-ahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yyn == 0 || yyn == YYTABLE_NINF)
	goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  /* Shift the look-ahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the token being shifted unless it is eof.  */
  if (yychar != YYEOF)
    yychar = YYEMPTY;

  *++yyvsp = yylval;


  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  yystate = yyn;
  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 2:
#line 114 "parser.yy"
    {if (interactive) std::cout << ". " << std::flush;}
    break;

  case 3:
#line 115 "parser.yy"
    {if (interactive) std::cout << ". " << std::flush;}
    break;

  case 4:
#line 118 "parser.yy"
    {}
    break;

  case 5:
#line 119 "parser.yy"
    {}
    break;

  case 6:
#line 120 "parser.yy"
    {if(interactive) yyerrok; else exit(1);}
    break;

  case 20:
#line 138 "parser.yy"
    {
    FILE *file = fopen(((yyvsp[0].stringptr))->c_str(), "r");
    if (!file) {
	std::cerr << "Failed to open file " << *((yyvsp[0].stringptr)) << std::endl;
    }
    else {
	console->checkModel(file);
	fclose(file);
    }
    delete (yyvsp[0].stringptr);
}
    break;

  case 21:
#line 151 "parser.yy"
    {
  std::string rngname;
  readRData((yyvsp[-1].pvec), _data_table, rngname);
  if (rngname.size() != 0) {
     std::cerr << "WARNING: .RNG.name assignment ignored" << std::endl;
  }
  delete (yyvsp[-1].pvec);
}
    break;

  case 22:
#line 159 "parser.yy"
    {}
    break;

  case 23:
#line 162 "parser.yy"
    {
    doDump(*(yyvsp[0].stringptr), DUMP_DATA, 1);
    delete (yyvsp[0].stringptr);
}
    break;

  case 24:
#line 168 "parser.yy"
    {
  if(open_data_buffer((yyvsp[0].stringptr))) {
    std::cout << "Reading data file " << *(yyvsp[0].stringptr) << std::endl;
  }
  else {
    std::cerr << "Unable to open file " << *(yyvsp[0].stringptr) << std::endl << std::flush;
  }
  delete (yyvsp[0].stringptr);
}
    break;

  case 25:
#line 179 "parser.yy"
    {
  std::cout << "Clearing data table " << std::endl;
  _data_table.clear();
}
    break;

  case 26:
#line 186 "parser.yy"
    {
    std::map<std::string, SArray> parameter_table;
    std::string rngname;
    readRData((yyvsp[-1].pvec), parameter_table, rngname);
    delete (yyvsp[-1].pvec);
    /* We have to set the name first, because the state or seed
       might be embedded in the parameter_table */
    if (rngname.size() != 0) {
        console->setRNGname(rngname, 1);
    }
    console->setParameters(parameter_table, 1);
}
    break;

  case 27:
#line 199 "parser.yy"
    {
    std::map<std::string, SArray> parameter_table;
    std::string rngname;
    readRData((yyvsp[-6].pvec), parameter_table, rngname);
    delete (yyvsp[-6].pvec);
    /* We have to set the name first, because the state or seed
       might be embedded in the parameter_table */
    if (rngname.size() != 0) {
        console->setRNGname(rngname, (yyvsp[-1].intval));
    }
    console->setParameters(parameter_table, (yyvsp[-1].intval));
}
    break;

  case 28:
#line 211 "parser.yy"
    {}
    break;

  case 29:
#line 214 "parser.yy"
    {
    doDump(*(yyvsp[0].stringptr), DUMP_PARAMETERS, 1);
    delete (yyvsp[0].stringptr);
}
    break;

  case 30:
#line 218 "parser.yy"
    {
    doDump(*(yyvsp[-5].stringptr), DUMP_PARAMETERS, (yyvsp[-1].intval));
    delete (yyvsp[-5].stringptr);
}
    break;

  case 31:
#line 224 "parser.yy"
    {
  if(open_data_buffer((yyvsp[0].stringptr))) {
    std::cout << "Reading parameter file " << *(yyvsp[0].stringptr) << std::endl;
  }
  else {
    std::cerr << "Unable to open file " << *(yyvsp[0].stringptr) << std::endl << std::flush;
  }
  delete (yyvsp[0].stringptr);
}
    break;

  case 32:
#line 233 "parser.yy"
    {
  /* Legacy option to not break existing scripts */
  if(open_data_buffer((yyvsp[0].stringptr))) {
    std::cout << "Reading initial values file " << *(yyvsp[0].stringptr) << std::endl;
  }
  else {
    std::cerr << "Unable to open file " << *(yyvsp[0].stringptr) << std::endl << std::flush;
  }
  delete (yyvsp[0].stringptr);
}
    break;

  case 33:
#line 247 "parser.yy"
    {
  console->compile(_data_table, 1, true);
}
    break;

  case 34:
#line 250 "parser.yy"
    {
    console->compile(_data_table, (yyvsp[-1].intval), true);
}
    break;

  case 35:
#line 255 "parser.yy"
    {
  if (!console->initialize()) {
    errordump();
  }
}
    break;

  case 36:
#line 262 "parser.yy"
    {
  //fixme: put refresh option
  long refresh = (yyvsp[0].intval)/40;
  updatestar((yyvsp[0].intval), refresh, 40);
}
    break;

  case 37:
#line 267 "parser.yy"
    {
  updatestar((yyvsp[-5].intval),(yyvsp[-1].intval), 40);
}
    break;

  case 38:
#line 272 "parser.yy"
    {return 0;}
    break;

  case 39:
#line 275 "parser.yy"
    {
  (yyval.ptree) = new ParseTree(P_VAR); setName((yyval.ptree), (yyvsp[0].stringptr));
}
    break;

  case 40:
#line 278 "parser.yy"
    {
  (yyval.ptree) = new ParseTree(P_VAR); setName((yyval.ptree), (yyvsp[-3].stringptr));
  setParameters((yyval.ptree), (yyvsp[-1].pvec));
}
    break;

  case 41:
#line 284 "parser.yy"
    {
  (yyval.pvec) = new std::vector<ParseTree*>(1, (yyvsp[0].ptree)); 
}
    break;

  case 42:
#line 287 "parser.yy"
    {
  (yyval.pvec)=(yyvsp[-2].pvec); (yyval.pvec)->push_back((yyvsp[0].ptree));
}
    break;

  case 43:
#line 292 "parser.yy"
    {
  (yyval.ptree) = new ParseTree(P_RANGE); setParameters((yyval.ptree),(yyvsp[0].ptree));
}
    break;

  case 44:
#line 295 "parser.yy"
    {
  (yyval.ptree) = new ParseTree(P_RANGE); setParameters((yyval.ptree), (yyvsp[-2].ptree), (yyvsp[0].ptree));
}
    break;

  case 45:
#line 301 "parser.yy"
    {(yyval.ptree) = new ParseTree(P_VALUE); (yyval.ptree)->setValue((yyvsp[0].intval));}
    break;

  case 48:
#line 308 "parser.yy"
    { 
  setMonitor((yyvsp[0].ptree), 1); delete (yyvsp[0].ptree);
}
    break;

  case 49:
#line 311 "parser.yy"
    { 
  setMonitor((yyvsp[-5].ptree), (yyvsp[-1].intval)); delete (yyvsp[-5].ptree);
}
    break;

  case 50:
#line 314 "parser.yy"
    {
  setMonitor((yyvsp[0].ptree), 1); delete (yyvsp[0].ptree);
}
    break;

  case 51:
#line 317 "parser.yy"
    { 
  setMonitor((yyvsp[-5].ptree), (yyvsp[-1].intval)); delete (yyvsp[-5].ptree);
}
    break;

  case 52:
#line 322 "parser.yy"
    {
  clearMonitor((yyvsp[0].ptree)); delete (yyvsp[0].ptree);
}
    break;

  case 53:
#line 331 "parser.yy"
    { (yyval.stringptr) = (yyvsp[0].stringptr);}
    break;

  case 54:
#line 332 "parser.yy"
    { (yyval.stringptr) = (yyvsp[0].stringptr); }
    break;

  case 55:
#line 335 "parser.yy"
    {
  doCoda ((yyvsp[0].ptree), "CODA"); delete (yyvsp[0].ptree);
}
    break;

  case 56:
#line 338 "parser.yy"
    {
  doCoda ((yyvsp[-5].ptree), *(yyvsp[-1].stringptr)); delete (yyvsp[-5].ptree); delete (yyvsp[-1].stringptr);
}
    break;

  case 57:
#line 341 "parser.yy"
    {
  doAllCoda ("CODA"); 
}
    break;

  case 58:
#line 344 "parser.yy"
    {
  doAllCoda (*(yyvsp[-1].stringptr)); delete (yyvsp[-1].stringptr); 
}
    break;

  case 59:
#line 349 "parser.yy"
    { loadModule(*(yyvsp[0].stringptr)); }
    break;

  case 60:
#line 354 "parser.yy"
    {
  (yyval.pvec) = new std::vector<ParseTree*>(1, (yyvsp[0].ptree));
}
    break;

  case 61:
#line 357 "parser.yy"
    {
  (yyval.pvec) = (yyvsp[-1].pvec); (yyval.pvec)->push_back((yyvsp[0].ptree));
}
    break;

  case 62:
#line 360 "parser.yy"
    {
  (yyval.pvec) = (yyvsp[-2].pvec); (yyval.pvec)->push_back((yyvsp[0].ptree));
}
    break;

  case 63:
#line 365 "parser.yy"
    {
  (yyval.ptree) = (yyvsp[0].ptree); setName((yyval.ptree), (yyvsp[-2].stringptr));
}
    break;

  case 64:
#line 368 "parser.yy"
    {
  (yyval.ptree) = new ParseTree(P_ARRAY);
  setName((yyval.ptree), (yyvsp[-2].stringptr));
  setParameters((yyval.ptree), (yyvsp[0].ptree));
}
    break;

  case 65:
#line 373 "parser.yy"
    {
  /* Allow this for setting the NAME of the random number generator */
  (yyval.ptree) = new ParseTree(P_VAR); setName((yyval.ptree), (yyvsp[-2].stringptr));
  ParseTree *p = new ParseTree(P_VAR); setName(p, (yyvsp[0].stringptr));
  setParameters((yyval.ptree), p);
}
    break;

  case 66:
#line 381 "parser.yy"
    {
  (yyval.ptree) = new ParseTree(P_ARRAY); 
  if ((yyvsp[-1].ptree)) 
    setParameters((yyval.ptree), (yyvsp[-3].ptree), (yyvsp[-1].ptree));
  else
    setParameters((yyval.ptree), (yyvsp[-3].ptree));
}
    break;

  case 68:
#line 393 "parser.yy"
    {(yyval.ptree)=0;}
    break;

  case 70:
#line 395 "parser.yy"
    {(yyval.ptree)=(yyvsp[0].ptree);}
    break;

  case 71:
#line 398 "parser.yy"
    {
  (yyval.ptree) = (yyvsp[0].ptree);
}
    break;

  case 74:
#line 407 "parser.yy"
    {(yyval.ptree) = (yyvsp[-1].ptree);}
    break;

  case 75:
#line 410 "parser.yy"
    { 
  (yyval.ptree) = new ParseTree(P_VECTOR); 
  setParameters((yyval.ptree), (yyvsp[0].ptree));
}
    break;

  case 76:
#line 414 "parser.yy"
    {
  (yyval.ptree) = new ParseTree(P_VECTOR);
  setParameters((yyval.ptree), (yyvsp[-1].pvec));
}
    break;

  case 77:
#line 420 "parser.yy"
    {(yyval.pvec) = new std::vector<ParseTree*>(1, (yyvsp[0].ptree)); }
    break;

  case 78:
#line 421 "parser.yy"
    {(yyval.pvec) = (yyvsp[-2].pvec); (yyval.pvec)->push_back((yyvsp[0].ptree));}
    break;

  case 79:
#line 424 "parser.yy"
    {(yyval.ptree) = new ParseTree(P_VALUE); (yyval.ptree)->setValue((yyvsp[0].val));}
    break;

  case 80:
#line 425 "parser.yy"
    {(yyval.ptree) = new ParseTree(P_VALUE); (yyval.ptree)->setValue(JAGS_NA);}
    break;

  case 81:
#line 431 "parser.yy"
    {;}
    break;

  case 82:
#line 434 "parser.yy"
    {;}
    break;

  case 83:
#line 435 "parser.yy"
    {;}
    break;

  case 84:
#line 438 "parser.yy"
    {;}
    break;

  case 85:
#line 439 "parser.yy"
    {;}
    break;

  case 86:
#line 442 "parser.yy"
    {;}
    break;

  case 87:
#line 443 "parser.yy"
    {;}
    break;

  case 88:
#line 444 "parser.yy"
    {;}
    break;

  case 89:
#line 445 "parser.yy"
    {;}
    break;

  case 90:
#line 446 "parser.yy"
    {;}
    break;

  case 91:
#line 447 "parser.yy"
    {;}
    break;

  case 92:
#line 450 "parser.yy"
    {;}
    break;

  case 94:
#line 454 "parser.yy"
    {;}
    break;

  case 95:
#line 455 "parser.yy"
    {;}
    break;

  case 96:
#line 458 "parser.yy"
    {;}
    break;

  case 97:
#line 459 "parser.yy"
    {;}
    break;

  case 98:
#line 462 "parser.yy"
    {;}
    break;

  case 99:
#line 463 "parser.yy"
    {;}
    break;


      default: break;
    }

/* Line 1126 of yacc.c.  */
#line 1945 "parser.cc"

  yyvsp -= yylen;
  yyssp -= yylen;


  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;


  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if YYERROR_VERBOSE
      yyn = yypact[yystate];

      if (YYPACT_NINF < yyn && yyn < YYLAST)
	{
	  int yytype = YYTRANSLATE (yychar);
	  YYSIZE_T yysize0 = yytnamerr (0, yytname[yytype]);
	  YYSIZE_T yysize = yysize0;
	  YYSIZE_T yysize1;
	  int yysize_overflow = 0;
	  char *yymsg = 0;
#	  define YYERROR_VERBOSE_ARGS_MAXIMUM 5
	  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
	  int yyx;

#if 0
	  /* This is so xgettext sees the translatable formats that are
	     constructed on the fly.  */
	  YY_("syntax error, unexpected %s");
	  YY_("syntax error, unexpected %s, expecting %s");
	  YY_("syntax error, unexpected %s, expecting %s or %s");
	  YY_("syntax error, unexpected %s, expecting %s or %s or %s");
	  YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s");
#endif
	  char *yyfmt;
	  char const *yyf;
	  static char const yyunexpected[] = "syntax error, unexpected %s";
	  static char const yyexpecting[] = ", expecting %s";
	  static char const yyor[] = " or %s";
	  char yyformat[sizeof yyunexpected
			+ sizeof yyexpecting - 1
			+ ((YYERROR_VERBOSE_ARGS_MAXIMUM - 2)
			   * (sizeof yyor - 1))];
	  char const *yyprefix = yyexpecting;

	  /* Start YYX at -YYN if negative to avoid negative indexes in
	     YYCHECK.  */
	  int yyxbegin = yyn < 0 ? -yyn : 0;

	  /* Stay within bounds of both yycheck and yytname.  */
	  int yychecklim = YYLAST - yyn;
	  int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
	  int yycount = 1;

	  yyarg[0] = yytname[yytype];
	  yyfmt = yystpcpy (yyformat, yyunexpected);

	  for (yyx = yyxbegin; yyx < yyxend; ++yyx)
	    if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
	      {
		if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
		  {
		    yycount = 1;
		    yysize = yysize0;
		    yyformat[sizeof yyunexpected - 1] = '\0';
		    break;
		  }
		yyarg[yycount++] = yytname[yyx];
		yysize1 = yysize + yytnamerr (0, yytname[yyx]);
		yysize_overflow |= yysize1 < yysize;
		yysize = yysize1;
		yyfmt = yystpcpy (yyfmt, yyprefix);
		yyprefix = yyor;
	      }

	  yyf = YY_(yyformat);
	  yysize1 = yysize + yystrlen (yyf);
	  yysize_overflow |= yysize1 < yysize;
	  yysize = yysize1;

	  if (!yysize_overflow && yysize <= YYSTACK_ALLOC_MAXIMUM)
	    yymsg = (char *) YYSTACK_ALLOC (yysize);
	  if (yymsg)
	    {
	      /* Avoid sprintf, as that infringes on the user's name space.
		 Don't have undefined behavior even if the translation
		 produced a string with the wrong number of "%s"s.  */
	      char *yyp = yymsg;
	      int yyi = 0;
	      while ((*yyp = *yyf))
		{
		  if (*yyp == '%' && yyf[1] == 's' && yyi < yycount)
		    {
		      yyp += yytnamerr (yyp, yyarg[yyi++]);
		      yyf += 2;
		    }
		  else
		    {
		      yyp++;
		      yyf++;
		    }
		}
	      yyerror (yymsg);
	      YYSTACK_FREE (yymsg);
	    }
	  else
	    {
	      yyerror (YY_("syntax error"));
	      goto yyexhaustedlab;
	    }
	}
      else
#endif /* YYERROR_VERBOSE */
	yyerror (YY_("syntax error"));
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse look-ahead token after an
	 error, discard it.  */

      if (yychar <= YYEOF)
        {
	  /* Return failure if at end of input.  */
	  if (yychar == YYEOF)
	    YYABORT;
        }
      else
	{
	  yydestruct ("Error: discarding", yytoken, &yylval);
	  yychar = YYEMPTY;
	}
    }

  /* Else will try to reuse look-ahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (0)
     goto yyerrorlab;

yyvsp -= yylen;
  yyssp -= yylen;
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (yyn != YYPACT_NINF)
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
	    {
	      yyn = yytable[yyn];
	      if (0 < yyn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
	YYABORT;


      yydestruct ("Error: popping", yystos[yystate], yyvsp);
      YYPOPSTACK;
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  *++yyvsp = yylval;


  /* Shift the error token. */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#ifndef yyoverflow
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEOF && yychar != YYEMPTY)
     yydestruct ("Cleanup: discarding lookahead",
		 yytoken, &yylval);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
		  yystos[*yyssp], yyvsp);
      YYPOPSTACK;
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
  return yyresult;
}


#line 466 "parser.yy"


int zzerror (const char *s)
{
  return_to_main_buffer();
  std::cerr << s << std::endl;
  return 0;
}

static Range getRange(ParseTree const *var)
{
  /* 
     Blank arguments, e.g. foo[] or bar[,1]  are not allowed.
  */
  unsigned int size = var->parameters().size();

  Index ind_lower(size), ind_upper(size);
  for (unsigned int i = 0; i < size; ++i) {
    ParseTree const *range_element = var->parameters()[i];
    switch(range_element->parameters().size()) {
    case 1:
      ind_lower[i] = (long) (range_element->parameters()[0]->value() + 1.0E-6);
      ind_upper[i] = ind_lower[i];
      break;
    case 2:
      ind_lower[i] = (long) (range_element->parameters()[0]->value() + 1.0E-6);  
      ind_upper[i] = (long) (range_element->parameters()[1]->value() + 1.0E-6);
      break;
    default:
      //Error! FIXME
      break;
    }
  }
  return Range(ind_lower, ind_upper);
}

void setMonitor(ParseTree const *var, int thin)
{
    std::string const &name = var->name();
    if (var->parameters().empty()) {
	/* Requesting the whole node */
	console->setMonitor(name, Range(), thin);
    }
    else {
	/* Requesting subset of a multivariate node */
	console->setMonitor(name, getRange(var), thin);
    }
}

void clearMonitor(ParseTree const *var)
{
    std::string const &name = var->name();
    if (var->parameters().empty()) {
	/* Requesting the whole node */
	console->clearMonitor(name, Range());
    }
    else {
	/* Requesting subset of a multivariate node */
	console->clearMonitor(name, getRange(var));
    }
}

/*
void doCoda (ParseTree const *var, std::string const &stem)
{
  unsigned int nchain = console->nchain();

  for (unsigned int n = 0; n < nchain; ++n) {
      std::ostringstream outstream;
    outstream << stem << n << ".out";
    std::string outname = outstream.str();
    std::ofstream out(outname.c_str());
    if (!out) {
      std::cerr << "Failed to open file " << outname << std::endl;
      return;
    }
    std::ostringstream indstream;
    indstream << stem << n << ".ind";
    std::string indname = indstream.str();
    std::ofstream ind(indname.c_str());
    if (!ind) {
      std::cerr << "Failed to open file " << indname << std::endl;
      out.close();
      return;
    }

    if (var->parameters().empty()) {
      console->coda(var->name(), Range(), out, ind, n);
    }
    else {
      console->coda(var->name(), getRange(var), out, ind, n);
    }
    out.close();
    ind.close();
  }
}
*/

void doAllCoda (std::string const &stem)
{

    /* Open index file */
    std::ostringstream indstream;
    indstream << stem << "index.txt";
    std::string indname = indstream.str();
    std::ofstream ind(indname.c_str());
    if (!ind) {
      std::cerr << "Failed to open file " << indname << std::endl;
      return;
    }

    std::vector<std::ofstream*> output;
    for (unsigned int n = 0; n < console->nchain(); ++n) {
	/* Open output file */
	std::ostringstream outstream;
	outstream << stem << "chain" << n + 1 << ".txt";
	std::string outname = outstream.str();
	std::ofstream *out = new std::ofstream(outname.c_str());
        if (*out) {
            output.push_back(out);
        }
	else {
	    ind.close();
	    for (unsigned int i = 0; i < output.size(); ++i) {
		output[i]->close();
	    }
	    std::cerr << "Failed to open file " << outname << std::endl;
	    return;
	}
    }
    console->coda(ind, output);
    ind.close();
    for (unsigned int i = 0; i < output.size(); ++i) {
	output[i]->close();
        delete output[i];
    }
}

void doCoda (ParseTree const *var, std::string const &stem)
{

    /* Open index file */
    std::ostringstream indstream;
    indstream << stem << "index.txt";
    std::string indname = indstream.str();
    std::ofstream ind(indname.c_str());
    if (!ind) {
      std::cerr << "Failed to open file " << indname << std::endl;
      return;
    }

    std::vector<std::ofstream*> output;
    for (unsigned int n = 0; n < console->nchain(); ++n) {
	/* Open output file */
	std::ostringstream outstream;
	outstream << stem << "chain" << n + 1 << ".txt";
	std::string outname = outstream.str();
	std::ofstream *out = new std::ofstream(outname.c_str());

        if (*out) {
            output.push_back(out);
        }
        else {
            ind.close();
            for (unsigned int i = 0; i < output.size(); ++i) {
                output[i]->close();
            }
            std::cerr << "Failed to open file " << outname << std::endl;
            return;
        }
    }

    if (var->parameters().empty()) {
      /* Requesting the whole node */
	console->coda(var->name(), Range(), ind, output);
    }
    else {
      /* Requesting subset of a multivariate node */
	console->coda(var->name(), getRange(var), ind, output);
    }

    ind.close();
    for (unsigned int i = 0; i < output.size(); ++i) {
	output[i]->close();
        delete output[i];
    }
}

/* Helper function for doDump that handles all the special cases
   (missing values etc) when writing a double value */
static void writeValue(double x, std::ostream &out, bool isdiscrete)
{
  using namespace std;

  if (x == JAGS_NA) {
    out << "NA";
  }
  else if (isnan(x)) {
    out << "NaN";
  }
  else if (isinf(x)) {
    if (x > 0)
      out << "Inf";
    else
      out << "-Inf";
  }
  else if (isdiscrete) {
      out << static_cast<int>(x);
  }
  else {
    out << x;
  }
}

void doDump(std::string const &file, DumpType type, unsigned int chain)
{
  std::map<std::string,SArray> data_table;
  std::string rng_name;
  if (!console->dumpState(data_table, rng_name, type, chain)) {
    return;
  }

  /* Open output file */
  std::ofstream out(file.c_str());
  if (!out) {
    std::cerr << "Failed to open file " << file << std::endl;
    return;
  }
  
  if (rng_name.size() != 0) {
    out << "\".RNG.name\" <- \"" << rng_name << "\"\n";
  }

  for (std::map<std::string, SArray>::const_iterator p = data_table.begin();
       p != data_table.end(); ++p) {
    std::string const &name = p->first;
    SArray const &sarray = p->second;
    double const *value = sarray.value();
    long length = sarray.length();
    out << "\"" << name << "\" <- " << std::endl;
    Index const &dim = sarray.dim(false);
    bool discrete = sarray.isDiscreteValued();

    if (discrete) {
       out << "as.integer(";
    }
    if (dim.size() == 1) {
      // Vector 
      if (dim[0] == 1) {
	// Scalar
	writeValue(value[0], out, sarray.isDiscreteValued());
      }
      else {
	// Vector of length > 1
	out << "c(";
	for (int i = 0; i < length; ++i) {
	  if (i > 0) {
	    out << ",";
	  }
	  writeValue(value[i], out, sarray.isDiscreteValued());
	}
	out << ")";
      }
    }
    else {
      // Array 
      out << "structure(c(";
      for (int i = 0; i < length; ++i) {
	if (i > 0) {
	  out << ",";
	}
	writeValue(value[i], out, sarray.isDiscreteValued());
      }
      out << "), .Dim = as.integer(c(";
      for (unsigned int j = 0; j < dim.size(); ++j) {
	if (j > 0) {
	  out << ",";
	}
	out << dim[j];
      }
      out << ")))";
    }
    if (discrete) 
      out << ")";
    out << "\n";
  }
  out.close();
}  

void setParameters(ParseTree *p, std::vector<ParseTree*> *parameters)
{
  /* 
     The parser dynamically allocates vectors of (pointers to)
     parameters. These vectors must be deleted when we are done with
     them.
  */
  p->setParameters(*parameters);
  delete parameters; 
}

void setParameters(ParseTree *p, ParseTree *param1)
{
  /*
    Wrapper function that creates a vector containing param1
    to be passed to ParseTree::setParameters.
  */
  std::vector<ParseTree *> parameters(1, param1);
  p->setParameters(parameters);
}

void setParameters(ParseTree *p, ParseTree *param1, ParseTree *param2)
{
  /*
    Wrapper function that creates a vector containing param1
    and param2, to be passed to ParseTree::setParameters
  */
  std::vector<ParseTree *> parameters;
  parameters.push_back(param1);
  parameters.push_back(param2);
  p->setParameters(parameters);
}

void setName(ParseTree *p, std::string *name)
{
  p->setName(*name);
  delete name;
}

static void errordump()
{
  if (console->model()) {
    std::ostringstream fname;
    for (unsigned int i = 1; i <= console->nchain(); ++i) {
       fname << "jags.dump" << i << ".R";
       std::cout << "Dumping chain " << i << " at iteration " 
                 << console->iter(i) << " to file " << fname.str() << std::endl;
       doDump(fname.str(), DUMP_ALL, i);
       fname.str("");
    }
  }
}

static void updatestar(long niter, long refresh, int width)
{
  if (refresh == 0) {
    console->update(niter);
    return;
  }

  if (width > niter / refresh + 1)
    width = niter / refresh + 1;

  std::cout << "Updating " << niter << std::endl;
  for (int i = 0; i < width - 1; ++i) {
    std::cout << "-";
  }
  std::cout << "| " << std::min(width * refresh, niter) << std::endl 
	    << std::flush;

  int col = 0;
  for (long n = niter; n > 0; n -= refresh) {
    long nupdate = std::min(n, refresh);
    if(console->update(nupdate))
      std::cout << "*" << std::flush;
    else {
      std::cout << std::endl;
      errordump();
      return;
    }
    col++;
    if (col == width || n <= nupdate) {
      int percent = 100 - (n-nupdate) * 100/niter;
      std::cout << " " << percent << "%" << std::endl;
      if (n > nupdate) {
	col = 0;
      }
    }
  }
}

static void loadModule(std::string const &name)
{
  std::cout << "   " << name << std::endl;
  lt_dlhandle mod = lt_dlopenext(name.c_str());
  if (mod == NULL) {
      std::cout << lt_dlerror() << std::endl;
  }
}

int main (int argc, char **argv)
{
  extern FILE *zzin;

  FILE *cmdfile = 0;
  if (argc > 2) {
    std::cerr << "Too many arguments" << std::endl;
  }
  else if (argc == 2) {
    interactive = false;
    cmdfile = fopen(argv[1],"r");
    if (cmdfile) {
      zzin = cmdfile;
    }
    else {
      std::cerr << "Unable to open command file " << argv[1] << std::endl;
      return 1;
    }
  }

  if(lt_dlinit()) {
      std::cerr << lt_dlerror() << std::endl;
      return 1;
  }

  /*
  pt2Func load_base = (pt2Func)(lt_dlsym(base, "load"));
  if (load_base == NULL) {
      std::cout << lt_dlerror() << std::endl;
      return 1;
  }
  else{
      (*load_base)();
  }
  */
  
  time_t t;
  time(&t);
  std::cout << "Welcome to " << PACKAGE_STRING << " on " << ctime(&t);
  std::cout << "JAGS is free software and comes with ABSOLUTELY NO WARRANTY" 
            << std::endl;

  std::cout << "Loading modules:" << std::endl;
  if (getenv("JAGS_HOME")) {
     /* Locate configuration file containing names of modules to load */
     std::string conf_name = std::string(getenv("JAGS_HOME")) + 
                             "/etc/modules.conf";
     std::ifstream conf_file(conf_name.c_str());
     if (conf_file) {
        std::string module_name;
        while((conf_file>> module_name)) {
           loadModule(module_name);
        }
     }
  }
  else {
    std::cerr << "Cannot find JAGS_HOME" << std::endl;
  }

  console = new Console(std::cout, std::cerr);

  zzparse();
  if (argc==2) {
    fclose(cmdfile);
  }

  delete console;
  lt_dlexit();
}

