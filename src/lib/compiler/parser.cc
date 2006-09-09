/* A Bison parser, made by GNU Bison 2.0.  */

/* Skeleton parser for Yacc-like parsing with Bison,
   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004 Free Software Foundation, Inc.

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
   Foundation, Inc., 59 Temple Place - Suite 330,
   Boston, MA 02111-1307, USA.  */

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

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Using locations.  */
#define YYLSP_NEEDED 0



/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     VAR = 258,
     DATA = 259,
     MODEL = 260,
     NAME = 261,
     FUNC = 262,
     IN = 263,
     ARROW = 264,
     FOR = 265,
     DOUBLE = 266,
     INT = 267,
     NEG = 268,
     ENDL = 269
   };
#endif
#define VAR 258
#define DATA 259
#define MODEL 260
#define NAME 261
#define FUNC 262
#define IN 263
#define ARROW 264
#define FOR 265
#define DOUBLE 266
#define INT 267
#define NEG 268
#define ENDL 269




/* Copy the first part of user declarations.  */
#line 2 "parser.yy"

#include <compiler/ParseTree.h>
#include "remap.h"
#include "parser_extra.h"

#include <cstdio>
#include <cfloat>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>
  
  void yyerror(const char *);
  int yylex();
  
#define YYDEBUG 1
#define YYERROR_VERBOSE 1
  
  static std::vector<ParseTree *> * _pvariables = 0;
  static ParseTree *_pdata = 0;
  static ParseTree *_prelations = 0;
  
  static void setName(ParseTree *p, std::string *name);
  static void setParameters(ParseTree *p, std::vector<ParseTree*> *parameters);
  static void setParameters(ParseTree *p, ParseTree *param1);
  static void setParameters(ParseTree *p, ParseTree *param1, 
			    ParseTree *param2);
  static void setParameters(ParseTree *p, ParseTree *param1, ParseTree *param2,
			    ParseTree *param3);
  static ParseTree *Truncated (ParseTree *left, ParseTree *right);
  

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

#if ! defined (YYSTYPE) && ! defined (YYSTYPE_IS_DECLARED)
#line 36 "parser.yy"
typedef union YYSTYPE {
  /*
    We can't put any C++ class that requires a constructor in
    union. In particular, we can't have strings.  To get round this
    limitation, we put a pointer to a string in the union, and then
    get the lexical analyzer to return a pointer to a dynamically
    allocated string.  The parser frees the memory when it is done.
    This is tedious, but the alternative is to copy character
    arrays using, e.g. strncpy instead.
  */
  double val;
  std::string *stringptr;
  ParseTree *ptree;
  std::vector<ParseTree*> *pvec;
} YYSTYPE;
/* Line 190 of yacc.c.  */
#line 152 "parser.cc"
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 213 of yacc.c.  */
#line 164 "parser.cc"

#if ! defined (yyoverflow) || YYERROR_VERBOSE

# ifndef YYFREE
#  define YYFREE free
# endif
# ifndef YYMALLOC
#  define YYMALLOC malloc
# endif

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   else
#    define YYSTACK_ALLOC alloca
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning. */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
# else
#  if defined (__STDC__) || defined (__cplusplus)
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   define YYSIZE_T size_t
#  endif
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
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
	  register YYSIZE_T yyi;		\
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
#define YYFINAL  13
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   134

/* YYNTOKENS -- Number of terminals. */
#define YYNTOKENS  30
/* YYNNTS -- Number of nonterminals. */
#define YYNNTS  25
/* YYNRULES -- Number of rules. */
#define YYNRULES  60
/* YYNRULES -- Number of states. */
#define YYNSTATES  118

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   269

#define YYTRANSLATE(YYX) 						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const unsigned char yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      25,    26,    15,    13,    20,    14,     2,    16,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    28,    19,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,    29,     2,     2,     2,     2,     2,
       2,    21,     2,    22,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    23,     2,    24,    27,     2,     2,     2,
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
       5,     6,     7,     8,     9,    10,    11,    12,    17,    18
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const unsigned char yyprhs[] =
{
       0,     0,     3,     4,     6,     9,    12,    16,    19,    23,
      25,    29,    31,    36,    38,    42,    47,    52,    56,    58,
      61,    63,    65,    67,    70,    73,    76,    83,    87,    94,
      98,   103,   107,   111,   115,   119,   121,   123,   128,   130,
     134,   136,   140,   143,   147,   149,   153,   155,   159,   160,
     162,   166,   171,   178,   184,   190,   195,   197,   199,   202,
     204
};

/* YYRHS -- A `-1'-separated list of the rules' RHS. */
static const yysigned_char yyrhs[] =
{
      31,     0,    -1,    -1,    37,    -1,    32,    37,    -1,    36,
      37,    -1,    32,    36,    37,    -1,     3,    33,    -1,     3,
      33,    19,    -1,    34,    -1,    33,    20,    34,    -1,     6,
      -1,     6,    21,    35,    22,    -1,    47,    -1,    35,    20,
      47,    -1,     4,    23,    39,    24,    -1,     5,    23,    39,
      24,    -1,    23,    39,    24,    -1,    40,    -1,    39,    40,
      -1,    44,    -1,    43,    -1,    41,    -1,    44,    19,    -1,
      43,    19,    -1,    42,    38,    -1,    10,    25,     6,     8,
      50,    26,    -1,    54,     9,    47,    -1,     7,    25,    54,
      26,     9,    47,    -1,    54,    27,    51,    -1,    54,    27,
      51,    52,    -1,    47,    15,    47,    -1,    45,    15,    47,
      -1,    47,    13,    47,    -1,    46,    13,    47,    -1,    54,
      -1,    11,    -1,     7,    25,    48,    26,    -1,    45,    -1,
      47,    16,    47,    -1,    46,    -1,    47,    14,    47,    -1,
      14,    47,    -1,    25,    47,    26,    -1,    47,    -1,    48,
      20,    47,    -1,    50,    -1,    49,    20,    50,    -1,    -1,
      47,    -1,    47,    28,    47,    -1,     7,    25,    48,    26,
      -1,    29,    25,    53,    20,    53,    26,    -1,    29,    25,
      20,    53,    26,    -1,    29,    25,    53,    20,    26,    -1,
      29,    25,    20,    26,    -1,    54,    -1,    11,    -1,    14,
      11,    -1,     6,    -1,     6,    21,    49,    22,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const unsigned short int yyrline[] =
{
       0,    83,    83,    84,    85,    86,    87,    90,    91,    94,
      95,    98,   101,   107,   108,   111,   118,   125,   131,   132,
     135,   136,   137,   138,   139,   142,   148,   154,   158,   174,
     178,   184,   189,   199,   204,   213,   214,   215,   219,   223,
     227,   231,   235,   239,   242,   243,   246,   247,   250,   251,
     252,   257,   264,   265,   266,   267,   270,   271,   272,   275,
     278
};
#endif

#if YYDEBUG || YYERROR_VERBOSE
/* YYTNME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals. */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "VAR", "DATA", "MODEL", "NAME", "FUNC",
  "IN", "ARROW", "FOR", "DOUBLE", "INT", "'+'", "'-'", "'*'", "'/'", "NEG",
  "ENDL", "';'", "','", "'['", "']'", "'{'", "'}'", "'('", "')'", "'~'",
  "':'", "'T'", "$accept", "input", "var_stmt", "dec_list", "node_dec",
  "dim_list", "data_stmt", "model_stmt", "relations", "relation_list",
  "relation", "for_loop", "counter", "determ_relation", "stoch_relation",
  "product", "sum", "expression", "expression_list", "range_list",
  "range_element", "distribution", "truncated", "param", "var", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const unsigned short int yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,    43,    45,    42,    47,   268,   269,    59,
      44,    91,    93,   123,   125,    40,    41,   126,    58,    84
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const unsigned char yyr1[] =
{
       0,    30,    31,    31,    31,    31,    31,    32,    32,    33,
      33,    34,    34,    35,    35,    36,    37,    38,    39,    39,
      40,    40,    40,    40,    40,    41,    42,    43,    43,    44,
      44,    45,    45,    46,    46,    47,    47,    47,    47,    47,
      47,    47,    47,    47,    48,    48,    49,    49,    50,    50,
      50,    51,    52,    52,    52,    52,    53,    53,    53,    54,
      54
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const unsigned char yyr2[] =
{
       0,     2,     0,     1,     2,     2,     3,     2,     3,     1,
       3,     1,     4,     1,     3,     4,     4,     3,     1,     2,
       1,     1,     1,     2,     2,     2,     6,     3,     6,     3,
       4,     3,     3,     3,     3,     1,     1,     4,     1,     3,
       1,     3,     2,     3,     1,     3,     1,     3,     0,     1,
       3,     4,     6,     5,     5,     4,     1,     1,     2,     1,
       4
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const unsigned char yydefact[] =
{
       2,     0,     0,     0,     0,     0,     0,     3,    11,     7,
       9,     0,     0,     1,     0,     4,     5,     0,     8,     0,
      59,     0,     0,     0,    18,    22,     0,    21,    20,     0,
       0,     6,     0,    36,     0,     0,     0,    38,    40,    13,
      35,    10,    48,     0,     0,    15,    19,     0,    25,    24,
      23,     0,     0,    16,     0,    42,     0,     0,    12,     0,
       0,     0,     0,     0,     0,    49,     0,    46,     0,     0,
       0,    27,     0,    29,    44,     0,    43,    14,    32,    34,
      33,    41,    31,    39,     0,    48,    60,     0,    48,    17,
       0,     0,    30,     0,    37,    50,    47,     0,     0,     0,
       0,    45,    28,    26,    51,    57,     0,     0,     0,    56,
      58,    55,     0,     0,    53,    54,     0,    52
};

/* YYDEFGOTO[NTERM-NUM]. */
static const yysigned_char yydefgoto[] =
{
      -1,     4,     5,     9,    10,    36,     6,     7,    48,    23,
      24,    25,    26,    27,    28,    37,    38,    65,    75,    66,
      67,    73,    92,   108,    40
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -47
static const yysigned_char yypact[] =
{
      94,    -1,     6,     8,    17,   103,    50,   -47,    37,    90,
     -47,    80,    80,   -47,    50,   -47,   -47,    19,   -47,    -1,
      60,    63,    75,    28,   -47,   -47,    61,    95,    97,     0,
      47,   -47,    88,   -47,    19,    19,    81,   100,   104,    59,
     -47,   -47,    19,   112,   113,   -47,   -47,    80,   -47,   -47,
     -47,    19,   114,   -47,    19,   -47,    79,    19,   -47,    19,
      19,    19,    19,    19,    19,    52,    84,   -47,    98,   115,
      54,    59,   101,    91,    59,    -2,   -47,    59,   -47,    96,
      96,    96,   -47,   -47,    19,    19,   -47,   116,    19,   -47,
      19,   102,   -47,    19,   -47,    59,   -47,    19,   105,    36,
      71,    59,    59,   -47,   -47,   -47,   111,    -4,   108,   -47,
     -47,   -47,   106,     2,   -47,   -47,   107,   -47
};

/* YYPGOTO[NTERM-NUM].  */
static const yysigned_char yypgoto[] =
{
     -47,   -47,   -47,   -47,   110,   -47,   125,     9,   -47,    -6,
     -19,   -47,   -47,   -47,   -47,   -47,   -47,   -14,    44,   -47,
     -46,   -47,   -47,   -44,   -11
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -1
static const unsigned char yytable[] =
{
      29,    29,    20,    39,    46,     8,    30,   105,    20,    51,
     106,    46,    29,   105,    15,    16,   106,    13,    93,    29,
      55,    56,   111,    31,    94,    20,    32,    52,   115,    11,
      33,    12,    68,    34,    20,    21,    29,    71,    22,    96,
      74,    70,    98,    77,    35,    78,    79,    80,    81,    82,
      83,    46,    45,    20,    21,     3,    93,    22,    17,    29,
      20,    21,   104,   112,    22,    61,    62,    63,    64,   116,
      95,    53,    61,    62,    63,    64,    74,    20,    89,   101,
      84,    42,   105,   102,    47,   106,    20,    21,    43,   109,
      22,   107,    61,    62,    63,    64,   109,     1,     2,     3,
      44,    57,   109,    58,    85,    76,    86,     2,     3,    18,
      19,    63,    64,    54,    49,    59,    50,    60,    20,    69,
      91,    72,   110,    88,    87,    97,    90,   100,   113,    41,
      14,   103,   114,   117,    99
};

static const unsigned char yycheck[] =
{
      11,    12,     6,    17,    23,     6,    12,    11,     6,     9,
      14,    30,    23,    11,     5,     6,    14,     0,    20,    30,
      34,    35,    26,    14,    26,     6,     7,    27,    26,    23,
      11,    23,    43,    14,     6,     7,    47,    51,    10,    85,
      54,    47,    88,    57,    25,    59,    60,    61,    62,    63,
      64,    70,    24,     6,     7,     5,    20,    10,    21,    70,
       6,     7,    26,   107,    10,    13,    14,    15,    16,   113,
      84,    24,    13,    14,    15,    16,    90,     6,    24,    93,
      28,    21,    11,    97,    23,    14,     6,     7,    25,   100,
      10,    20,    13,    14,    15,    16,   107,     3,     4,     5,
      25,    20,   113,    22,    20,    26,    22,     4,     5,    19,
      20,    15,    16,    25,    19,    15,    19,    13,     6,     6,
      29,     7,    11,     8,    26,     9,    25,    25,    20,    19,
       5,    26,    26,    26,    90
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const unsigned char yystos[] =
{
       0,     3,     4,     5,    31,    32,    36,    37,     6,    33,
      34,    23,    23,     0,    36,    37,    37,    21,    19,    20,
       6,     7,    10,    39,    40,    41,    42,    43,    44,    54,
      39,    37,     7,    11,    14,    25,    35,    45,    46,    47,
      54,    34,    21,    25,    25,    24,    40,    23,    38,    19,
      19,     9,    27,    24,    25,    47,    47,    20,    22,    15,
      13,    13,    14,    15,    16,    47,    49,    50,    54,     6,
      39,    47,     7,    51,    47,    48,    26,    47,    47,    47,
      47,    47,    47,    47,    28,    20,    22,    26,     8,    24,
      25,    29,    52,    20,    26,    47,    50,     9,    50,    48,
      25,    47,    47,    26,    26,    11,    14,    20,    53,    54,
      11,    26,    53,    20,    26,    26,    53,    26
};

#if ! defined (YYSIZE_T) && defined (__SIZE_TYPE__)
# define YYSIZE_T __SIZE_TYPE__
#endif
#if ! defined (YYSIZE_T) && defined (size_t)
# define YYSIZE_T size_t
#endif
#if ! defined (YYSIZE_T)
# if defined (__STDC__) || defined (__cplusplus)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# endif
#endif
#if ! defined (YYSIZE_T)
# define YYSIZE_T unsigned int
#endif

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
    { 								\
      yyerror ("syntax error: cannot back up");\
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
      yysymprint (stderr, 					\
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
  unsigned int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %u), ",
             yyrule - 1, yylno);
  /* Print the symbols being reduced, and their result.  */
  for (yyi = yyprhs[yyrule]; 0 <= yyrhs[yyi]; yyi++)
    YYFPRINTF (stderr, "%s ", yytname [yyrhs[yyi]]);
  YYFPRINTF (stderr, "-> %s\n", yytname [yyr1[yyrule]]);
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
   SIZE_MAX < YYSTACK_BYTES (YYMAXDEPTH)
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
  register const char *yys = yystr;

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
  register char *yyd = yydest;
  register const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

#endif /* !YYERROR_VERBOSE */



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
  
  register int yystate;
  register int yyn;
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
  register short int *yyssp;

  /* The semantic value stack.  */
  YYSTYPE yyvsa[YYINITDEPTH];
  YYSTYPE *yyvs = yyvsa;
  register YYSTYPE *yyvsp;



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


  yyvsp[0] = yylval;

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
	yyoverflow ("parser stack overflow",
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),

		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyoverflowlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyoverflowlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	short int *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyoverflowlab;
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
        case 7:
#line 90 "parser.yy"
    { _pvariables = (yyvsp[0].pvec); }
    break;

  case 8:
#line 91 "parser.yy"
    { _pvariables = (yyvsp[-1].pvec); }
    break;

  case 9:
#line 94 "parser.yy"
    { (yyval.pvec) = new std::vector<ParseTree*>(1, (yyvsp[0].ptree)); }
    break;

  case 10:
#line 95 "parser.yy"
    { (yyval.pvec) = (yyvsp[-2].pvec); (yyval.pvec)->push_back((yyvsp[0].ptree)); }
    break;

  case 11:
#line 98 "parser.yy"
    {
  (yyval.ptree) = new ParseTree(P_VAR); setName((yyval.ptree), (yyvsp[0].stringptr)); 
}
    break;

  case 12:
#line 101 "parser.yy"
    {
  (yyval.ptree) = new ParseTree(P_VAR); setName((yyval.ptree), (yyvsp[-3].stringptr));
  setParameters((yyval.ptree), (yyvsp[-1].pvec));
}
    break;

  case 13:
#line 107 "parser.yy"
    { (yyval.pvec) = new std::vector<ParseTree*>(1, (yyvsp[0].ptree)); }
    break;

  case 14:
#line 108 "parser.yy"
    { (yyval.pvec)=(yyvsp[-2].pvec); (yyval.pvec)->push_back((yyvsp[0].ptree)); }
    break;

  case 15:
#line 111 "parser.yy"
    {
  ParseTree *p = new ParseTree(P_RELATIONS);
  setParameters(p, (yyvsp[-1].pvec));
  _pdata = p;
}
    break;

  case 16:
#line 118 "parser.yy"
    { 
  ParseTree *p = new ParseTree(P_RELATIONS);
  setParameters(p, (yyvsp[-1].pvec));
  _prelations = p;
}
    break;

  case 17:
#line 125 "parser.yy"
    {
  (yyval.ptree) = new ParseTree(P_RELATIONS);
  setParameters((yyval.ptree), (yyvsp[-1].pvec));
}
    break;

  case 18:
#line 131 "parser.yy"
    { (yyval.pvec) = new std::vector<ParseTree*>(1, (yyvsp[0].ptree)); }
    break;

  case 19:
#line 132 "parser.yy"
    { (yyval.pvec)=(yyvsp[-1].pvec); (yyval.pvec)->push_back((yyvsp[0].ptree)); }
    break;

  case 25:
#line 142 "parser.yy"
    {
  (yyval.ptree) = new ParseTree(P_FOR);
  setParameters((yyval.ptree), (yyvsp[-1].ptree), (yyvsp[0].ptree));
}
    break;

  case 26:
#line 148 "parser.yy"
    {
  (yyval.ptree) = new ParseTree(P_COUNTER); setName((yyval.ptree), (yyvsp[-3].stringptr));
  setParameters((yyval.ptree), (yyvsp[-1].ptree));
}
    break;

  case 27:
#line 154 "parser.yy"
    {
  (yyval.ptree) = new ParseTree(P_DETRMREL);
  setParameters((yyval.ptree), (yyvsp[-2].ptree), (yyvsp[0].ptree));
}
    break;

  case 28:
#line 158 "parser.yy"
    {

  /* 
     The link function is given using an S-style replacement function
     notation.  We need to turn this round so the inverse link
     function is applied to the RHS of the deterministic relation
  */
     
  ParseTree *p = new ParseTree(P_LINK); setName(p, (yyvsp[-5].stringptr));
  setParameters(p, (yyvsp[0].ptree));

  (yyval.ptree) = new ParseTree(P_DETRMREL);
  setParameters((yyval.ptree), (yyvsp[-3].ptree), p);
}
    break;

  case 29:
#line 174 "parser.yy"
    {
  (yyval.ptree) = new ParseTree(P_STOCHREL); 
  setParameters((yyval.ptree), (yyvsp[-2].ptree), (yyvsp[0].ptree));
}
    break;

  case 30:
#line 178 "parser.yy"
    {
  (yyval.ptree) = new ParseTree(P_STOCHREL); 
  setParameters((yyval.ptree), (yyvsp[-3].ptree), (yyvsp[-1].ptree), (yyvsp[0].ptree));
}
    break;

  case 31:
#line 184 "parser.yy"
    {
  (yyval.pvec) = new std::vector<ParseTree*>;
  (yyval.pvec)->push_back((yyvsp[-2].ptree));
  (yyval.pvec)->push_back((yyvsp[0].ptree));
}
    break;

  case 32:
#line 189 "parser.yy"
    {
  //  This creates a shift-reduce conflict because in the expression
  //  A*B*C, (A*B) is a valid expression. By default, bison shifts,
  //  which is what we want. The warning is suppressed with the %expect
  //  declaration (See also sum: below).
  (yyval.pvec) = (yyvsp[-2].pvec);
  (yyval.pvec)->push_back((yyvsp[0].ptree));
}
    break;

  case 33:
#line 199 "parser.yy"
    {
  (yyval.pvec) = new std::vector<ParseTree*>;
  (yyval.pvec)->push_back((yyvsp[-2].ptree));
  (yyval.pvec)->push_back((yyvsp[0].ptree));
}
    break;

  case 34:
#line 204 "parser.yy"
    {
  //  This creates a shift-reduce conflict. By default, bison shifts,
  //  which is what we want. The warning is suppressed with the %expect
  //  declaration (See also product: above).
  (yyval.pvec) = (yyvsp[-2].pvec);
  (yyval.pvec)->push_back((yyvsp[0].ptree));
}
    break;

  case 36:
#line 214 "parser.yy"
    {(yyval.ptree) = new ParseTree(P_VALUE); (yyval.ptree)->setValue((yyvsp[0].val));}
    break;

  case 37:
#line 215 "parser.yy"
    {
  (yyval.ptree) = new ParseTree(P_FUNCTION); setName((yyval.ptree), (yyvsp[-3].stringptr));
  setParameters((yyval.ptree), (yyvsp[-1].pvec));
}
    break;

  case 38:
#line 219 "parser.yy"
    {
  (yyval.ptree) = new ParseTree(P_OPERATOR); (yyval.ptree)->setOperator(OP_MULTIPLY);
  setParameters((yyval.ptree), (yyvsp[0].pvec));
}
    break;

  case 39:
#line 223 "parser.yy"
    {
  (yyval.ptree) = new ParseTree(P_OPERATOR); (yyval.ptree)->setOperator(OP_DIVIDE);
  setParameters((yyval.ptree), (yyvsp[-2].ptree), (yyvsp[0].ptree));
}
    break;

  case 40:
#line 227 "parser.yy"
    {
  (yyval.ptree) = new ParseTree(P_OPERATOR); (yyval.ptree)->setOperator(OP_ADD);
  setParameters((yyval.ptree), (yyvsp[0].pvec));
}
    break;

  case 41:
#line 231 "parser.yy"
    {
  (yyval.ptree) = new ParseTree(P_OPERATOR); (yyval.ptree)->setOperator(OP_SUBTRACT);
  setParameters((yyval.ptree), (yyvsp[-2].ptree), (yyvsp[0].ptree));
}
    break;

  case 42:
#line 235 "parser.yy"
    {
  (yyval.ptree) = new ParseTree(P_OPERATOR); (yyval.ptree)->setOperator(OP_NEG);
  setParameters((yyval.ptree), (yyvsp[0].ptree));
}
    break;

  case 43:
#line 239 "parser.yy"
    { (yyval.ptree) = (yyvsp[-1].ptree); }
    break;

  case 44:
#line 242 "parser.yy"
    { (yyval.pvec) = new std::vector<ParseTree*>(1, (yyvsp[0].ptree)); }
    break;

  case 45:
#line 243 "parser.yy"
    { (yyval.pvec)=(yyvsp[-2].pvec); (yyval.pvec)->push_back((yyvsp[0].ptree)); }
    break;

  case 46:
#line 246 "parser.yy"
    { (yyval.pvec) = new std::vector<ParseTree*>(1, (yyvsp[0].ptree)); }
    break;

  case 47:
#line 247 "parser.yy"
    { (yyval.pvec)=(yyvsp[-2].pvec); (yyval.pvec)->push_back((yyvsp[0].ptree)); }
    break;

  case 48:
#line 250 "parser.yy"
    {(yyval.ptree) = new ParseTree(P_RANGE);}
    break;

  case 49:
#line 251 "parser.yy"
    {(yyval.ptree) = new ParseTree(P_RANGE); setParameters((yyval.ptree),(yyvsp[0].ptree));}
    break;

  case 50:
#line 252 "parser.yy"
    {
  (yyval.ptree) = new ParseTree(P_RANGE); setParameters((yyval.ptree), (yyvsp[-2].ptree), (yyvsp[0].ptree));
}
    break;

  case 51:
#line 258 "parser.yy"
    {
  (yyval.ptree) = new ParseTree(P_DENSITY); setName((yyval.ptree), (yyvsp[-3].stringptr));
  setParameters((yyval.ptree), (yyvsp[-1].pvec));
}
    break;

  case 52:
#line 264 "parser.yy"
    {(yyval.ptree) = Truncated((yyvsp[-3].ptree),(yyvsp[-1].ptree));}
    break;

  case 53:
#line 265 "parser.yy"
    {(yyval.ptree) = Truncated(0,(yyvsp[-1].ptree));}
    break;

  case 54:
#line 266 "parser.yy"
    {(yyval.ptree) = Truncated((yyvsp[-2].ptree),0);}
    break;

  case 55:
#line 267 "parser.yy"
    {(yyval.ptree) = Truncated(0,0);}
    break;

  case 57:
#line 271 "parser.yy"
    {(yyval.ptree) = new ParseTree(P_VALUE); (yyval.ptree)->setValue((yyvsp[0].val));}
    break;

  case 58:
#line 272 "parser.yy"
    {(yyval.ptree) = new ParseTree(P_VALUE); (yyval.ptree)->setValue(-(yyvsp[0].val));}
    break;

  case 59:
#line 275 "parser.yy"
    {
  (yyval.ptree) = new ParseTree(P_VAR); setName((yyval.ptree), (yyvsp[0].stringptr));
}
    break;

  case 60:
#line 278 "parser.yy"
    {
  (yyval.ptree) = new ParseTree(P_VAR); setName((yyval.ptree), (yyvsp[-3].stringptr));
  setParameters((yyval.ptree), (yyvsp[-1].pvec));
}
    break;


    }

/* Line 1037 of yacc.c.  */
#line 1502 "parser.cc"

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
	  YYSIZE_T yysize = 0;
	  int yytype = YYTRANSLATE (yychar);
	  const char* yyprefix;
	  char *yymsg;
	  int yyx;

	  /* Start YYX at -YYN if negative to avoid negative indexes in
	     YYCHECK.  */
	  int yyxbegin = yyn < 0 ? -yyn : 0;

	  /* Stay within bounds of both yycheck and yytname.  */
	  int yychecklim = YYLAST - yyn;
	  int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
	  int yycount = 0;

	  yyprefix = ", expecting ";
	  for (yyx = yyxbegin; yyx < yyxend; ++yyx)
	    if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
	      {
		yysize += yystrlen (yyprefix) + yystrlen (yytname [yyx]);
		yycount += 1;
		if (yycount == 5)
		  {
		    yysize = 0;
		    break;
		  }
	      }
	  yysize += (sizeof ("syntax error, unexpected ")
		     + yystrlen (yytname[yytype]));
	  yymsg = (char *) YYSTACK_ALLOC (yysize);
	  if (yymsg != 0)
	    {
	      char *yyp = yystpcpy (yymsg, "syntax error, unexpected ");
	      yyp = yystpcpy (yyp, yytname[yytype]);

	      if (yycount < 5)
		{
		  yyprefix = ", expecting ";
		  for (yyx = yyxbegin; yyx < yyxend; ++yyx)
		    if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
		      {
			yyp = yystpcpy (yyp, yyprefix);
			yyp = yystpcpy (yyp, yytname[yyx]);
			yyprefix = " or ";
		      }
		}
	      yyerror (yymsg);
	      YYSTACK_FREE (yymsg);
	    }
	  else
	    yyerror ("syntax error; also virtual memory exhausted");
	}
      else
#endif /* YYERROR_VERBOSE */
	yyerror ("syntax error");
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse look-ahead token after an
	 error, discard it.  */

      if (yychar <= YYEOF)
        {
          /* If at end of input, pop the error token,
	     then the rest of the stack, then return failure.  */
	  if (yychar == YYEOF)
	     for (;;)
	       {

		 YYPOPSTACK;
		 if (yyssp == yyss)
		   YYABORT;
		 yydestruct ("Error: popping",
                             yystos[*yyssp], yyvsp);
	       }
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

#ifdef __GNUC__
  /* Pacify GCC when the user code never invokes YYERROR and the label
     yyerrorlab therefore never appears in user code.  */
  if (0)
     goto yyerrorlab;
#endif

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
  yydestruct ("Error: discarding lookahead",
              yytoken, &yylval);
  yychar = YYEMPTY;
  yyresult = 1;
  goto yyreturn;

#ifndef yyoverflow
/*----------------------------------------------.
| yyoverflowlab -- parser overflow comes here.  |
`----------------------------------------------*/
yyoverflowlab:
  yyerror ("parser stack overflow");
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
  return yyresult;
}


#line 284 "parser.yy"


using std::cout;
using std::endl;

void yyerror (const char *s)
{
  cout << s << endl;
}

static ParseTree *Truncated (ParseTree *left, ParseTree *right)
{
  ParseTree *p = new ParseTree(P_BOUNDS);
  setParameters(p, left, right);
  return p;
}

void setName(ParseTree *p, std::string *name)
{
  /* 
     The scanner cannot return a string, because a string cannot be
     part of a union.  So instead the scanner returns a pointer to
     a newly allocated string.  This memory must be freed as soon
     as it is used.
  */
  p->setName(*name);
  delete name;
}
	

void setParameters(ParseTree *p, std::vector<ParseTree*> *parameters)
{
  /* 
     Same as setName (above).  The parser dynamically allocates
     vectors of (pointers to) parameters. These vectors must be
     deleted when we are done with them.
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

void setParameters(ParseTree *p, ParseTree *param1, ParseTree *param2,
		   ParseTree *param3)
{
  /*
    Wrapper function that creates a vector containing param1, param2
    and param3, to be passed to ParseTree::setParameters
  */
  std::vector<ParseTree *> parameters;
  parameters.push_back(param1);
  parameters.push_back(param2);
  parameters.push_back(param3);
  p->setParameters(parameters);
}

int parse_bugs (FILE *file, std::vector<ParseTree*> **dec_list, 
                ParseTree **data, ParseTree **relations)
{
  extern FILE *yyin;
  yyin = file;

  delete _pvariables; _pvariables = 0;
  delete _prelations; _prelations = 0;
  delete _pdata; _pdata = 0;
  if (yyparse() == 0) {
      *dec_list = _pvariables;
      *data = _pdata;
      *relations = _prelations;
      return 0;
  }
  else {
    // FIXME: memory leaks when we have a syntax error?
    // FIXME TOO: yylineno is not remapped. Remove it
    extern int yylineno;
    return yylineno;
  }
}

