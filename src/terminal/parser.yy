/* -*-C++-*- */
%{
#include <config.h>

#ifdef Win32
#include <windows.h>   /* For getCurrentDirectory */
#include <io.h>        /* For chdir */
#else
#include <unistd.h>    /* For getcwd, chdir */
#endif

//#include <limits.h>

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <fstream>
#include <list>
#include <iterator>
#include <string>
#include <utility>

#include <dirent.h>
#include <time.h>
#include <errno.h>

#include <Console.h>
#include <module/Module.h>
#include <compiler/ParseTree.h>
#include <util/nainf.h>
#include <cstring>
#include <ltdl.h>

//Required for warning about masked distributions after module loading
#include <deque>
#include <distribution/Distribution.h>
#include <compiler/Compiler.h>

#include "ReadData.h"

    typedef void(*pt2Func)();

    int zzerror(const char *);
    int zzlex();
    int zzlex_destroy();
#define YYERROR_VERBOSE 0
    static jags::Console *console;
    bool interactive;
    extern int command_buffer_count;
    void setName(jags::ParseTree *p, std::string *name);
    std::map<std::string, jags::SArray> _data_table;
    std::deque<lt_dlhandle> _dyn_lib;
    bool open_data_buffer(std::string const *name);
    bool open_command_buffer(std::string const *name);
    void return_to_main_buffer();
    void setMonitor(jags::ParseTree const *var, int thin, std::string const &type);
    void clearMonitor(jags::ParseTree const *var, std::string const &type);
    void doCoda (jags::ParseTree const *var, std::string const &stem, std::string const &type);
    void doAllCoda (std::string const &stem, std::string const &type);
    void dumpNodeNames (std::string const &file, std::string const &type);
    void doDump (std::string const &file, jags::DumpType type, unsigned int chain);
    void dumpMonitors(std::string const &file, std::string const &type);
    void doSystem(std::string const *command);
    std::string ExpandFileName(char const *s);

    static bool getWorkingDirectory(std::string &name);
    static void errordump();
    static void updatestar(long niter, long refresh, int width);
	// Run adaptation phase until adapted, regardless of iterations:
    static void autoadaptstar(long maxiter, long refresh, int width);
    static void adaptstar(long niter, long refresh, int width, bool force);
    static void setParameters(jags::ParseTree *p, jags::ParseTree *param1);
    static void setParameters(jags::ParseTree *p, std::vector<jags::ParseTree*> *parameters);
    static void setParameters(jags::ParseTree *p, jags::ParseTree *param1, jags::ParseTree *param2);
    static void loadModule(std::string const &name);
    static void unloadModule(std::string const &name);
    static void dumpSamplers(std::string const &file);
    static void delete_pvec(std::vector<jags::ParseTree*> *);
    static void print_unused_variables(std::map<std::string, jags::SArray> const &table, bool data);
    static void listFactories(jags::FactoryType type);
	static void listModules();
    static void setFactory(std::string const &name, jags::FactoryType type,
                           std::string const &status);
    static void setSeed(unsigned int seed);
    static bool Jtry(bool ok);
	// Needed for update (and adapt) functions to dump variable states:
    static bool Jtry_dump(bool ok);
	// Note: on !ok, Jtry(_dump) only exits if !interactive and otherwise returns ok, 
	// so the parent function should evaluate the value and return; if appriopriate
    %}

%defines
%name-prefix "zz"
%union
{
  int intval;
  double val;
  std::string *stringptr;
  jags::ParseTree *ptree;
  std::vector<jags::ParseTree*> *pvec;
  std::vector<double> *vec;
  std::vector<long> *ivec;
}

%token <intval> INT
%token <val>    DOUBLE
%token <stringptr> NAME
%token <stringptr> STRING
%token <stringptr> SYSCMD
%token <intval> ENDCMD

%token <intval> MODEL
%token <intval> DATA
%token <intval> IN 
%token <intval> TO
%token <intval> INITS
%token <intval> PARAMETERS

%token <intval> COMPILE
%token <intval> INITIALIZE
%token <intval> ADAPT
%token <intval> AUTOADAPT
%token <intval> FORCEADAPT
%token <intval> UPDATE
%token <intval> BY
%token <intval> MONITORS
%token <intval> NODENAMES
%token <intval> MONITOR
%token <intval> TYPE
%token <intval> SET
%token <intval> CLEAR
%token <intval> THIN
%token <intval> CODA
%token <intval> STEM
%token <intval> EXIT
%token <intval> NCHAINS
%token <intval> CHAIN
%token <intval> LOAD
%token <intval> UNLOAD
%token <intval> SAMPLER
%token <intval> SAMPLERS
%token <intval> RNGTOK
%token <intval> FACTORY;
%token <intval> FACTORIES;
%token <intval> MODULES;
%token <intval> SEED;

%token <intval> LIST 
%token <intval> STRUCTURE
%token <intval> DIM
%token <intval> NA
%token <intval> R_NULL
%token <intval> DIMNAMES
%token <intval> ITER
%token <intval> ARROW
%token <intval> ENDDATA
%token <intval> ASINTEGER
%token <intval> DOTDATA

%token <intval> DIRECTORY
%token <intval> CD
%token <intval> PWD
%token <intval> RUN
%token <intval> ENDSCRIPT

%type <ptree> var index 
%type <ptree> r_assignment r_structure
%type <ptree> range_element r_dim
%type <ptree> r_attribute_list
%type <ptree> r_value 
%type <pvec>  r_value_list r_assignment_list range_list
%type <ptree> r_value_collection r_integer_collection r_collection r_data
%type <stringptr> file_name;
%type <stringptr> r_name;

%%

input: {
    if (interactive && command_buffer_count == 0) 
	std::cout << ". " << std::flush;
}
| input line {
    if (interactive && command_buffer_count == 0) 
	std::cout << ". " << std::flush;
}
;

line: ENDCMD {}
| command ENDCMD {}
| error ENDCMD {if(interactive) yyerrok; else exit(1); }
| run_script {}
| ENDSCRIPT ENDCMD {}
| SYSCMD ENDCMD { doSystem($1); delete $1;}
;

command: model 
| data_in
| data_to
| data_clear
| parameters_in
| parameters_to
| compile
| initialize
| adapt
| autoadapt
| forceadapt
| update
| monitor
| monitors_to
| node_names_to
| coda
| load
| unload
| exit
| read_dir
| get_working_dir
| set_working_dir
| samplers_to
| list_factories
| list_modules
| set_factory
| set_seed
;

model: MODEL IN file_name {
    std::FILE *file = std::fopen(ExpandFileName(($3)->c_str()).c_str(), "r");
    if (!file) {
	std::cerr << "Failed to open file " << *($3) << std::endl;
	if (!interactive) exit(1);
    }
    else {
	Jtry(console->checkModel(file));
	std::fclose(file);
    }
    delete $3;
 }
| MODEL CLEAR {
    console->clearModel();
 }
;

data_in: data r_assignment_list ENDDATA {
    std::string rngname;
    readRData($2, _data_table, rngname);
    if (rngname.size() != 0) {
	std::cerr << "WARNING: .RNG.name assignment ignored" << std::endl;
    }
    delete_pvec($2);
 }
| data {
    // Failed to open the data file 
    if (!interactive) exit(1);
  }    
;

data_to: DATA TO file_name {
    doDump(*$3, jags::DUMP_DATA, 1);
    delete $3;
}
;

data: DATA IN file_name {
    if(open_data_buffer($3)) {
	std::cout << "Reading data file " << *$3 << std::endl;
    }
    else {
	std::cerr << "Unable to open file " << *$3 << std::endl << std::flush;
	if (!interactive) exit(1);
    }
    delete $3;
 }
;

data_clear: DATA CLEAR {
    std::cout << "Clearing data table " << std::endl;
    _data_table.clear();
}
;

parameters_in: parameters r_assignment_list ENDDATA
{
    std::map<std::string, jags::SArray> parameter_table;
    std::string rngname;
    readRData($2, parameter_table, rngname);
    delete_pvec($2);
    /* Set all chains to the same state. If the user sets the
       RNG state in addition to the parameter values then all
       chains will be identical!
    */
    if (console->model() == 0) {
	std::cout << "ERROR: Initial values ignored. "
		  <<  "(You must compile the model first)" << std::endl;
	if (!interactive) exit(1);
    }
    for (unsigned int i = 1; i <= console->nchain(); ++i) {
	/* We have to set the name first, because the state or seed
	   might be embedded in the parameter_table */
	if (rngname.size() != 0) {
	    Jtry(console->setRNGname(rngname, i));
	}
	Jtry(console->setParameters(parameter_table, i));
    }
    print_unused_variables(parameter_table, false);

}
| parameters r_assignment_list ENDDATA ',' CHAIN '(' INT ')' 
{
    std::map<std::string, jags::SArray> parameter_table;
    std::string rngname;
    readRData($2, parameter_table, rngname);
    delete $2;
    /* We have to set the name first, because the state or seed
       might be embedded in the parameter_table */
    if (rngname.size() != 0) {
        Jtry(console->setRNGname(rngname, $7));
    }
    Jtry(console->setParameters(parameter_table, $7));
    print_unused_variables(parameter_table, false);
}
| parameters {} // Failed to open the file
;

parameters_to: PARAMETERS TO file_name {
    doDump(*$3, jags::DUMP_PARAMETERS, 1);
    delete $3;
}
| PARAMETERS TO file_name ',' CHAIN '(' INT ')' {
    doDump(*$3, jags::DUMP_PARAMETERS, $7);
    delete $3;
}
;

parameters: PARAMETERS IN file_name {
  if(open_data_buffer($3)) {
    std::cout << "Reading parameter file " << *$3 << std::endl;
  }
  else {
    std::cerr << "Unable to open file " << *$3 << std::endl << std::flush;
  }
  delete $3;
}
| INITS IN file_name {
  /* Legacy option to not break existing scripts */
  if(open_data_buffer($3)) {
    std::cout << "Reading initial values file " << *$3 << std::endl;
  }
  else {
    std::cerr << "Unable to open file " << *$3 << std::endl << std::flush;
  }
  delete $3;
}
;

compile: COMPILE {
    Jtry(console->compile(_data_table, 1, true));
    print_unused_variables(_data_table, true);
 }
| COMPILE ',' NCHAINS '(' INT ')' {
    Jtry(console->compile(_data_table, $5, true));
    print_unused_variables(_data_table, true);
}
;

initialize: INITIALIZE {
    if (!console->initialize()) {
	errordump();
    }
}
;

autoadapt: AUTOADAPT INT {
    long refresh = interactive ? $2/50 : 0;
    autoadaptstar($2, refresh, 50);
}
| AUTOADAPT INT ',' BY '(' INT ')' {
    autoadaptstar($2,$6, 50);
}
;

forceadapt: FORCEADAPT INT {
    long refresh = interactive ? $2/50 : 0;
    adaptstar($2, refresh, 50, true);
}
| FORCEADAPT INT ',' BY '(' INT ')' {
    adaptstar($2,$6, 50, true);
}
;

adapt: ADAPT INT {
    long refresh = interactive ? $2/50 : 0;
    adaptstar($2, refresh, 50, false);
}
| ADAPT INT ',' BY '(' INT ')' {
    adaptstar($2,$6, 50, false);
}
;

update: UPDATE INT {
    long refresh = interactive ? $2/50 : 0;
    updatestar($2, refresh, 50);
}
| UPDATE INT ',' BY '(' INT ')' {
  updatestar($2,$6, 50);
}
;

exit: EXIT { return 0; }
;

var: NAME {
  $$ = new jags::ParseTree(jags::P_VAR); setName($$, $1);
}
| NAME '[' range_list ']' {
  $$ = new jags::ParseTree(jags::P_VAR); setName($$, $1);
  setParameters($$, $3);
}
;

range_list: range_element {
  $$ = new std::vector<jags::ParseTree*>(1, $1); 
}
| range_list ',' range_element {
  $$=$1; $$->push_back($3);
}
;

range_element: index {
  $$ = new jags::ParseTree(jags::P_RANGE); setParameters($$, $1);
}
| index ':' index {
  $$ = new jags::ParseTree(jags::P_RANGE); setParameters($$, $1, $3);
}
;

/* FIXME: Use integer value here */
index: INT {$$ = new jags::ParseTree(jags::P_VALUE); $$->setValue($1);}
;

monitor: monitor_set
| monitor_clear
;

monitor_set: MONITOR SET var  { 
    setMonitor($3, 1, "trace"); delete $3;
}
| MONITOR SET var ',' THIN '(' INT ')' { 
    setMonitor($3, $7, "trace"); delete $3;
}
| MONITOR var {
    setMonitor($2, 1, "trace"); delete $2;
}
| MONITOR var ',' THIN '(' INT ')' { 
    setMonitor($2, $6, "trace"); delete $2;
}
| MONITOR var ',' TYPE '(' NAME ')' {
    setMonitor($2, 1, *$6);
    delete $6;
}
| MONITOR var ',' TYPE '(' NAME ')' THIN '(' INT ')' {
    setMonitor($2, $10, *$6); 
    delete $6;
}
| MONITOR var ',' THIN '(' INT ')' TYPE '(' NAME ')' {
    setMonitor($2, $6, *$10); 
    delete $10;
}
;

monitor_clear: MONITOR CLEAR var {
    clearMonitor($3, "trace"); delete $3;
}
| MONITOR CLEAR var ',' TYPE '(' NAME ')' {
    clearMonitor($3, *$7);
    delete $7;
}
;

monitors_to:  MONITORS TO file_name 
{
    dumpMonitors(*$3, "trace");
    delete $3;
}
|
MONITORS TO file_name ',' TYPE '(' NAME ')' {
    dumpMonitors(*$3, *$7);
    delete $3;
    delete $7; 
}
;

node_names_to:  NODENAMES TO file_name 
{
    dumpNodeNames(*$3, "all");
    delete $3;
}
|
NODENAMES TO file_name ',' TYPE '(' NAME ')' {
    dumpNodeNames(*$3, *$7);
    delete $3;
    delete $7; 
}
;


/* 
   File names may optionally be enclosed in quotes, and this is required
   if the name includes spaces.
*/
file_name: NAME { $$ = $1;}
| STRING { $$ = $1; }
;

coda: CODA var {
  doCoda ($2, "CODA", "*"); delete $2;
}
| CODA var ',' STEM '(' file_name ')' {
  doCoda ($2, *$6, "*"); delete $2; delete $6;
}
| CODA var ',' STEM '(' file_name ')' TYPE '(' NAME ')' {
  doCoda ($2, *$6, *$10); delete $2; delete $6; delete $10;
}
| CODA '*' {
  doAllCoda ("CODA", "*"); 
}
| CODA '*' ',' STEM '(' file_name ')' {
  doAllCoda (*$6, "*"); delete $6; 
}
| CODA '*' ',' STEM '(' file_name ')' TYPE '(' NAME ')' {
  doAllCoda (*$6, *$10); delete $6; delete $10;
}
;

load: LOAD file_name { loadModule(*$2); }
;

unload: UNLOAD NAME { unloadModule(*$2); }
;

samplers_to: SAMPLERS TO file_name 
{
    dumpSamplers(*$3);
    delete $3;
}
;

list_factories: LIST FACTORIES ',' TYPE '(' SAMPLER ')'
{
    listFactories(jags::SAMPLER_FACTORY);
}
|
LIST FACTORIES ',' TYPE '(' RNGTOK ')'
{
    listFactories(jags::RNG_FACTORY);
}
|
LIST FACTORIES ',' TYPE '(' MONITOR ')'
{
    listFactories(jags::MONITOR_FACTORY);
}
;

list_modules: LIST MODULES
{
    listModules();
}
;

set_factory: SET FACTORY STRING NAME ',' TYPE '(' SAMPLER ')'
{
    setFactory(*$3, jags::SAMPLER_FACTORY, *$4);
    delete $3;
    delete $4;
}
|
SET FACTORY NAME NAME ',' TYPE '(' RNGTOK ')'
{
    setFactory(*$3, jags::RNG_FACTORY, *$4);
    delete $3;
    delete $4;
}
|
SET FACTORY NAME NAME ',' TYPE '(' MONITOR ')'
{
    setFactory(*$3, jags::MONITOR_FACTORY, *$4);
    delete $3;
    delete $4;
}
;

set_seed: SET SEED INT
{
    setSeed($3);
}
;

/* Rules for scanning dumped R datasets */

r_assignment_list: r_assignment {
  $$ = new std::vector<jags::ParseTree*>(1, $1);
}
| r_assignment_list r_assignment {
  $$ = $1; $$->push_back($2);
}
| r_assignment_list ';' r_assignment {
  $$ = $1; $$->push_back($3);
}
;

r_assignment: r_name ARROW r_structure {
  $$ = $3; setName($$, $1);
}
| r_name ARROW r_collection {
  $$ = new jags::ParseTree(jags::P_ARRAY);
  setName($$, $1);
  setParameters($$, $3);
}
| r_name ARROW STRING {
  /* Allow this for setting the NAME of the random number generator */
  $$ = new jags::ParseTree(jags::P_VAR); setName($$, $1);
  jags::ParseTree *p = new jags::ParseTree(jags::P_VAR); setName(p, $3);
  setParameters($$, p);
}
;

r_name: STRING
| NAME
| '`' NAME '`' {
    /* R >= 2.4.0 uses backticks for quoted names */
    $$ = $2;
}

r_data: r_collection
| DOTDATA '=' r_collection {
    $$ = $3;
}

r_structure: STRUCTURE '(' r_data ',' r_attribute_list ')' {
  $$ = new jags::ParseTree(jags::P_ARRAY); 
  if ($5) 
    setParameters($$, $3, $5);
  else
    setParameters($$, $3);
}
| STRUCTURE '(' r_data ')' {
    $$ = new jags::ParseTree(jags::P_ARRAY);
    setParameters($$, $3);
}
;

/* The only attribute we are interested in is .Dim. The rest are
   simply discarded - see below */
r_attribute_list: r_dim
| r_generic_attribute {$$=0;}
| r_attribute_list ',' r_generic_attribute
| r_attribute_list ',' r_dim {$$=$3;}
;

r_dim: DIM '=' r_collection {
  $$ = $3;
}
| DIM '=' range_element {
  $$ = $3;
}
;

r_collection: r_integer_collection
| r_value_collection
;

r_integer_collection: ASINTEGER '(' r_value_collection ')' {$$ = $3;}
;

r_value_collection: r_value { 
  $$ = new jags::ParseTree(jags::P_VECTOR); 
  setParameters($$, $1);
}
| 'c' '(' r_value_list ')' {
  $$ = new jags::ParseTree(jags::P_VECTOR);
  setParameters($$, $3);
}
;
 
r_value_list: r_value {$$ = new std::vector<jags::ParseTree*>(1, $1); }
| r_value_list ',' r_value {$$ = $1; $$->push_back($3);}
;

r_value: DOUBLE {$$ = new jags::ParseTree(jags::P_VALUE); $$->setValue($1);}
| NA {$$ = new jags::ParseTree(jags::P_VALUE); $$->setValue(JAGS_NA);}
;

/* Rules for parsing generic attributes.  We don't want to do anything
   with the results, just have the parser accept them */

r_generic_attribute: NAME '=' r_generic_vector {;}
;

r_generic_list: r_generic_list_element {;}
| r_generic_list ',' r_generic_list_element {;}
;

r_generic_list_element: r_generic_vector {;}
| NAME '=' r_generic_vector {;}
;

r_generic_vector: r_numeric_vector {;}
| ASINTEGER '(' r_numeric_vector ')' {;}
| r_character_vector {;}
| LIST '(' r_generic_list ')' {;}
| STRUCTURE '(' r_generic_list ')' {;}
| R_NULL {;}
;

r_numeric_vector: DOUBLE {;}
| 'c' '(' r_double_list ')'
;

r_double_list: DOUBLE {;}
| r_double_list ',' DOUBLE {;}
;

r_character_vector: STRING {;}
| 'c' '(' r_string_list ')' {;}
;

r_string_list: STRING {;}
| r_string_list ',' STRING {;}
;

/* Rules for interacting with the operating system */

get_working_dir: PWD
{
    std::string name;
    if (getWorkingDirectory(name)) {
	std::cout << name << std::endl;
    }
    else {
	std::cout << "ERROR: " << name << std::endl;
    }
}

set_working_dir: CD file_name
{
    if (chdir(($2)->c_str()) == -1) {
	std::cout << "ERROR: Cannot change working directory" << std::endl;
    }
    delete $2;
}

read_dir: DIRECTORY
{
    std::string name;
    if (!getWorkingDirectory(name)) {
	std::cerr << "ERROR: Unable to get working directory name\n"
		  << name << std::endl;
	return 0;
    }
	
    DIR *dir;
    struct dirent *de;
    if ((dir = opendir(name.c_str())) != 0) {
	while ((de = readdir(dir)) != 0) {
	    if (std::strcmp(de->d_name, ".") && std::strcmp(de->d_name, "..")) {
		std::cout << de->d_name << "\n";
	    }
	}
	closedir(dir);
    }
    else {
	std::cerr << "Unable to open working directory" << std::endl;
    }
}    

run_script: RUN file_name {
    if(open_command_buffer($2)) {
	std::cout << "Running script file " << *$2 << std::endl;
    }
    else {
	std::cerr << "Unable to open script file " << *$2 << std::endl;
    }
    delete $2;
 }
;

%%

int zzerror (const char *s)
{
    return_to_main_buffer();
    std::cerr << s << std::endl;
    return 0;
}

static jags::Range getRange(jags::ParseTree const *var)
{
  /* 
     Blank arguments, e.g. foo[] or bar[,1]  are not allowed.
  */
  unsigned int size = var->parameters().size();

  std::vector<unsigned long>  ind_lower(size), ind_upper(size);
  for (unsigned int i = 0; i < size; ++i) {
    jags::ParseTree const *range_element = var->parameters()[i];
    switch(range_element->parameters().size()) {
    case 1:
	ind_lower[i] = static_cast<int>(range_element->parameters()[0]->value());
      ind_upper[i] = ind_lower[i];
      break;
    case 2:
	ind_lower[i] = static_cast<int>(range_element->parameters()[0]->value());
	ind_upper[i] = static_cast<int>(range_element->parameters()[1]->value());
      break;
    default:
      //Error! FIXME
      break;
    }
  }
  return jags::SimpleRange(ind_lower, ind_upper);
}

void setMonitor(jags::ParseTree const *var, int thin, std::string const &type)
{
    std::string const &name = var->name();
    if (var->parameters().empty()) {
	/* Requesting the whole node */
	console->setMonitor(name, jags::Range(), thin, type);
    }
    else {
	/* Requesting subset of a multivariate node */
	console->setMonitor(name, getRange(var), thin, type);
    }
}

void clearMonitor(jags::ParseTree const *var, std::string const &type)
{
    std::string const &name = var->name();
    if (var->parameters().empty()) {
	/* Requesting the whole node */
	console->clearMonitor(name, jags::Range(), type);
    }
    else {
	/* Requesting subset of a multivariate node */
	console->clearMonitor(name, getRange(var), type);
    }
}

void doAllCoda (std::string const &stem, std::string const &type)
{
    console->coda(stem, type);
}

void doCoda (jags::ParseTree const *var, std::string const &stem, std::string const &type)
{
    //FIXME: Allow list of several nodes

    std::vector<std::pair<std::string, jags::Range> > dmp;
    if (var->parameters().empty()) {
	/* Requesting the whole node */
	dmp.push_back(std::pair<std::string,jags::Range>(var->name(), jags::Range()));
    }
    else {
	/* Requesting subset of a multivariate node */
	dmp.push_back(std::pair<std::string,jags::Range>(var->name(), getRange(var)));
    }
    console->coda(dmp, stem, type);
}

/* Helper function for doDump that handles all the special cases
   (missing values etc) when writing a double value */
static void writeValue(double x, std::ostream &out, bool isdiscrete)
{
  using namespace std;

  if (x == JAGS_NA) {
    out << "NA";
  }
  else if (jags_isnan(x)) {
    out << "NaN";
  }
  else if (!jags_finite(x)) {
    if (x > 0)
      out << "Inf";
    else
      out << "-Inf";
  }
  else if (isdiscrete) {
      out << static_cast<int>(x) << "L";
  }
  else {
    out << x;
  }
}

void doDump(std::string const &file, jags::DumpType type, unsigned int chain)
{
    std::map<std::string,jags::SArray> data_table;
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
	out << "`.RNG.name` <- \"" << rng_name << "\"\n";
    }

    for (std::map<std::string, jags::SArray>::const_iterator p = data_table.begin();
	 p != data_table.end(); ++p) {
	std::string const &name = p->first;
	jags::SArray const &sarray = p->second;
	std::vector<double> const &value = sarray.value();
	long length = sarray.length();
	out << "`" << name << "` <- " << std::endl;
	std::vector<unsigned long> const &dim = sarray.dim(false);
	bool discrete = sarray.isDiscreteValued();

	if (dim.size() == 1) {
	    // Vector 
	    if (dim[0] == 1) {
		// Scalar
		writeValue(value[0], out, discrete);
	    }
	    else {
		// Vector of length > 1
		out << "c(";
		for (int i = 0; i < length; ++i) {
		    if (i > 0) {
			out << ",";
		    }
		    writeValue(value[i], out, discrete);
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
		writeValue(value[i], out, discrete);
	    }
	    out << "), .Dim = c(";
	    for (unsigned int j = 0; j < dim.size(); ++j) {
		if (j > 0) {
		    out << ",";
		}
		out << dim[j] << "L";
	    }
	    out << "))";
	}
	out << "\n";
    }
    out.close();
}  

void dumpMonitors(std::string const &file, std::string const &type)
{
    std::map<std::string,jags::SArray> data_table;

    if (!console->dumpMonitors(data_table, type, false)) {
	return;
    }

    /* Open output file */
    std::ofstream out(file.c_str());
    if (!out) {
	std::cerr << "Failed to open file " << file << std::endl;
	return;
    }

    out << "`" << type << "` <-\nstructure(list(";

    std::map<std::string, jags::SArray>::const_iterator p;
    for (p = data_table.begin(); p != data_table.end(); ++p) {
	std::string const &name = p->first;
	jags::SArray const &sarray = p->second;
	std::vector<double> const &value = sarray.value();
	long length = sarray.length();

	if (p != data_table.begin()) {
	    out << ", \n";
	}
	out << "\"" << name << "\" = ";
	std::vector<unsigned long> const &dim = sarray.dim(false);
	bool discrete = sarray.isDiscreteValued();
	bool named = !sarray.dimNames().empty();

	if (dim.size() == 1 && !named) {
	    // Vector 
	    if (dim[0] == 1) {
		// Scalar
		writeValue(value[0], out, discrete);
	    }
	    else {
		// Vector of length > 1
		out << "c(";
		for (int i = 0; i < length; ++i) {
		    if (i > 0) {
			out << ",";
		    }
		    writeValue(value[i], out, discrete);
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
		writeValue(value[i], out, discrete);
	    }
	    out << "), .Dim = ";
	    if (named) {
		out << "structure(";
	    }
	    out << "c(";
	    for (unsigned int j = 0; j < dim.size(); ++j) {
		if (j > 0) {
		    out << ",";
		}
		out << dim[j] << "L";
	    }
	    out << ")";
	    if (named) {
		std::vector<std::string> const &dnames = sarray.dimNames();
		out << ", .Names = c(";
		for (unsigned int k = 0; k < dnames.size(); ++k) {
		    if (k > 0) {
			out << ",";
		    }
		    out << "\"" << dnames[k] << "\"";
		}
		out << "))";
	    }
	    out << ")";
	}
    }

    out << "), \n.Names = c(";
    for (p = data_table.begin(); p != data_table.end(); ++p) {
	if (p != data_table.begin()) {
	    out << ", ";
	}
	std::string const &name = p->first;
	out << "\"" << name << "\"";
    }
    out << "))";
    out.close();
}

void dumpNodeNames(std::string const &file, std::string const &type)
{
    std::vector<std::string> node_names;

    /* Open output file */
    std::ofstream out(file.c_str());
    if (!out) {
	std::cerr << "Failed to open file " << file << std::endl;
	return;
    }

	std::vector<std::string> alltypes;
	
	// Special rule for "all" types:
	if ( type == "all" ) {
		alltypes.push_back("constant");
		alltypes.push_back("deterministic");
		alltypes.push_back("stochastic");
		alltypes.push_back("fixed");
		alltypes.push_back("observations");
	}
	else {
		alltypes.push_back(type);
	}
	
	for ( unsigned int i=0; i < alltypes.size(); i++ ){
	
	    out << "`" << alltypes[i] << "` <-\n";
	
		console->dumpNodeNames(node_names, alltypes[i], true);

		// If no matching names just write an empty character vector:
	    if ( node_names.size() == 0 ) {
			out << "character(0)\n" << std::endl;
			out.close();
			return;
	    }
	
		out << "c(";
		unsigned int r = 0;
	    for (unsigned int i = 0;  i < node_names.size(); i++) {
			out << "\"" << node_names[i] << "\"";
			if( (i+1) < node_names.size() ) {
				out << ", ";
			}
			r++;
			if( r == 5 ) {
				out << "\n";
				r = 0;
			}
		}
		out << ")\n\n";
	
	}
	
	out.close();
}

void setParameters(jags::ParseTree *p, std::vector<jags::ParseTree*> *parameters)
{
  /* 
     The parser dynamically allocates vectors of (pointers to)
     parameters. These vectors must be deleted when we are done with
     them.
  */
  p->setParameters(*parameters);
  delete parameters; 
}

void setParameters(jags::ParseTree *p, jags::ParseTree *param1)
{
  /*
    Wrapper function that creates a vector containing param1
    to be passed to jags::ParseTree::setParameters.
  */
  std::vector<jags::ParseTree *> parameters(1, param1);
  p->setParameters(parameters);
}

void setParameters(jags::ParseTree *p, jags::ParseTree *param1, jags::ParseTree *param2)
{
  /*
    Wrapper function that creates a vector containing param1
    and param2, to be passed to jags::ParseTree::setParameters
  */
  std::vector<jags::ParseTree *> parameters;
  parameters.push_back(param1);
  parameters.push_back(param2);
  p->setParameters(parameters);
}

void setName(jags::ParseTree *p, std::string *name)
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
		      << console->iter() << " to file " << fname.str() 
		      << std::endl;
	    doDump(fname.str(), jags::DUMP_ALL, i);
	    fname.str("");
	}
	// Moved clearModel from Console.cc CATCH_ERRORS to here
	// to allow doDump to work as described in the manual:
	console->clearModel();
    }
    if (!interactive) exit(1);
}

static void updatestar(long niter, long refresh, int width)
{
    std::cout << "Updating " << niter << std::endl;

    bool adapt = console->isAdapting();
    if (adapt && console->iter() > 0) {
	//Turn off iteration immediately if we have some burn-in
	if (console->adaptOff()) {
	    adapt = false;
	}
	else {
	    std::cout << std::endl;
	    errordump();
	    return;
	}
    }

    if (refresh == 0) {
	if( !Jtry_dump(console->update(niter/2)) ) return;
	bool status = true;
	if (adapt) {
	    if (!console->checkAdaptation(status)) {
		errordump();
		return;
	    }
	    if (!console->adaptOff()) {
		errordump();
		return;
	    }
	}
	if( !Jtry_dump(console->update(niter - niter/2)) ) return;
	if (!status) {
	    std::cerr << "WARNING: Adaptation incomplete\n";
	}
	return;
    }

    if (width > niter / refresh + 1)
	width = niter / refresh + 1;

    for (int i = 0; i < width - 1; ++i) {
	std::cout << "-";
    }
    std::cout << "| " << std::min(width * refresh, niter) << std::endl 
	      << std::flush;

    int col = 0;
    bool status = true;
    for (long n = niter; n > 0; n -= refresh) {
	if (adapt && n <= niter/2) {
	    // Turn off adaptive mode half way through burnin
	    if (!console->checkAdaptation(status)) {
		std::cout << std::endl;
		errordump();
		return;
	    }
	    if (console->adaptOff()) {
		adapt = false;
	    }
	    else {
		std::cout << std::endl;
		errordump();
		return;
	    }
	}
	long nupdate = std::min(n, refresh);
	if(Jtry_dump(console->update(nupdate))) {
	    std::cout << "*" << std::flush;
	}
	else {
	    std::cout << std::endl;
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
    if (!status) {
	std::cerr << "WARNING: Adaptation incomplete\n";
    }
}

static void autoadaptstar(long maxiter, long refresh, int width)
{
    if (!console->isAdapting()) {
        std::cout << "Adaptation skipped: model is not in adaptive mode.\n";
		return;
    }
    std::cout << "Autoadapting up to " << maxiter << " iterations" << std::endl;
	
    if (width > maxiter / refresh + 1)
		width = maxiter / refresh + 1;
	
	long i = 0;
	bool status = false;
	
	if ( refresh == 0 ) {
		for (i = 0; i < maxiter; i++) {
			if (!console->checkAdaptation(status)) {
			    errordump();
			    return;
			}
			if(status)
				break;

			if( !Jtry_dump(console->update(1)) ) {
				std::cout << std::endl;
				return;
			}
		}
	}
	else {
	    for (int j = 0; j < width - 1; ++j) {
			std::cout << "-";
	    }
	    std::cout << "| " << std::min(width * refresh, maxiter) << std::endl 
		      << std::flush;

		int col = 0;
	    for (long n = maxiter; n > 0; n -= refresh) {
		    long nupdate = std::min(n, refresh);
			for (long j = 0; j < nupdate; j++) {
				if (!console->checkAdaptation(status)) {
				    errordump();
				    return;
				}
				if(status)
					break;

				if( !Jtry_dump(console->update(1)) ) {
					std::cout << std::endl;
					return;
				}
				i++;
			}
			if(status){
		    	std::cout << std::endl;
				break;
			}

	        std::cout << "+" << std::flush;
	    	col++;
	    	if (col == width || n <= nupdate) {
	    	    int percent = 100 - (n-nupdate) * 100/maxiter;
	    	    std::cout << " " << percent << "%" << std::endl;
	    	    if (n > nupdate) {
	        		col = 0;
		        }
		    }
	    }
	}
		
	if (!console->checkAdaptation(status)) {
	    errordump();
	    return;
	}
	
	if (!status) {
	    std::cerr << "Adaptation incomplete\n";
	}
	else {
		std::cout << "Adaptation completed in " << i << " iterations" << std::endl;
		if (!console->adaptOff()) {
			std::cout << std::endl;
			errordump();
			return;
	    }
	}
    return;
}

static void adaptstar(long niter, long refresh, int width, bool force)
{
    if (!console->isAdapting()) {
	    if( force ) {
			// Missing endl is deliberate - updatestar will write to the same line:
		    std::cout << "Model not in adaptive mode: ";
			updatestar(niter, refresh, width);			
		}
		else {
	        std::cerr << "Adaptation skipped: model is not in adaptive mode.\n";
		}
		return;
    }
    std::cout << "Adapting " << niter << std::endl;
    
    bool status = true;
    if (refresh == 0) {
	if( !Jtry_dump(console->update(niter)) ) return;
	if (!console->checkAdaptation(status)) {
	    errordump();
	    return;
	}
	if (!status) {
	    std::cerr << "Adaptation incomplete\n";
	    return;
	}
	else {
	    std::cerr << "Adaptation successful\n";
	    return;
	}
    }

    if (width > niter / refresh + 1)
	width = niter / refresh + 1;

    for (int i = 0; i < width - 1; ++i) {
	std::cout << "-";
    }
    std::cout << "| " << std::min(width * refresh, niter) << std::endl 
	      << std::flush;

    int col = 0;
    for (long n = niter; n > 0; n -= refresh) {
	long nupdate = std::min(n, refresh);
	if(Jtry_dump(console->update(nupdate)))
	    std::cout << "+" << std::flush;
	else {
	    std::cout << std::endl;
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
    if (!console->checkAdaptation(status)) {
	std::cout << std::endl;
	errordump();
	return;
    }
    if (!status) {
	std::cerr << "Adaptation incomplete.\n";
    }
    else {
	std::cerr << "Adaptation successful\n";
    }
}

static void loadModule(std::string const &name)
{
    std::cout << "Loading module: " << name;
    lt_dlhandle mod = lt_dlopenext(name.c_str());
    if (mod == NULL) {
	std::cout << ": " << lt_dlerror() << std::endl;
    }
    else {
	std::cout << ": ok" << std::endl;
	_dyn_lib.push_front(mod);
	jags::Console::loadModule(name);
    }
}

static void unloadModule(std::string const &name)
{
    std::cout << "Unloading module: " << name << std::endl;
    jags::Console::unloadModule(name);
}

void exiting() {
	std::cout << "Exiting JAGS" << std::endl;
}

int main (int argc, char **argv)
{
  std::atexit(exiting);

  extern std::FILE *zzin;

  std::FILE *cmdfile = 0;
  if (argc > 2) {
    std::cerr << "Too many arguments" << std::endl;
  }
  else if (argc == 2) {
    interactive = false;
    cmdfile = std::fopen(ExpandFileName(argv[1]).c_str(),"r");
    if (cmdfile) {
      zzin = cmdfile;
    }
    else {
      std::cerr << "Unable to open command file " << argv[1] << std::endl;
      return 1;
    }
  }
  else {
    interactive = true;
  }

#ifndef _WIN32
  /* 
     - Allows emulation of dynamic loading on platforms that do not
     support it by preloading modules. 
     - Causes build failures on mingw-w64 (as at 21 April 2010) so
     not used on Windows platform.
  */
  LTDL_SET_PRELOADED_SYMBOLS();
#endif

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
  loadModule("basemod");
  loadModule("bugs");

  console = new jags::Console(std::cout, std::cerr);

  zzparse();
  zzlex_destroy();

  if (argc==2) {
      std::fclose(cmdfile);
  }
  
  //Unload modules
  std::vector<std::string> mods = jags::Console::listModules();
  for (unsigned int i = 0; i < mods.size(); ++i) {
      jags::Console::unloadModule(mods[i]);
  }
  delete console;
  //Release dynamic libraries. 
  for (unsigned int i = 0; i < _dyn_lib.size(); ++i) {
      lt_dlclose(_dyn_lib[i]);
  }
  
  lt_dlexit();
}

static bool getWorkingDirectory(std::string &name)
{
    char buf[FILENAME_MAX];
#ifdef Win32
    if (getCurrentDirectory(FILENAME_MAX, buf)) {
	name = buf;
	return true;
    }
    else {
	name = "Error in getCurrentDirectory";
	return false;
    }
#else
    if (getcwd(buf, FILENAME_MAX)) {
	name = buf;
	return true;
    }
    else {
	switch(errno) {
	case EACCES:
	    name = "Access denied";
	    break;
	case ENOENT:
	    name = "Not found";
	    break;
	case ERANGE:
	    name = "Directory name too long";
	    break;
	default:
	    name = "Error in getcwd";
	    break;
	}
	return false;
    }
#endif
}

static void dumpSamplers(std::string const &file)
{
    std::ofstream out(file.c_str());
    if (!out) {
	std::cerr << "Failed to open file " << file << std::endl;
	return;
    }

    std::vector<std::vector<std::string> > sampler_list;
    console->dumpSamplers(sampler_list);
    for (unsigned int i = 0; i < sampler_list.size(); ++i) {
	for (unsigned int j = 1; j < sampler_list[i].size(); ++j) {
	    out << i + 1 << "\t" 
		<< sampler_list[i][0] << "\t" //First element is sampler name
		<< sampler_list[i][j] << "\n"; //Rest are node names
	}
    }

    out.close();
}

static void delete_pvec(std::vector<jags::ParseTree*> *pv)
{
    for (unsigned int i = 0; i < pv->size(); ++i) {
	delete (*pv)[i];
    }
    delete pv;
}

static void print_unused_variables(std::map<std::string, jags::SArray> const &table,
				   bool data)
{
    std::vector<std::string> supplied_vars;
    for (std::map<std::string, jags::SArray>::const_iterator p = table.begin();
	 p != table.end(); ++p)
    {
	supplied_vars.push_back(p->first);
    }
    
    std::vector<std::string> unused_vars;
    std::vector<std::string> model_vars = console->variableNames();
    if (!data) {
	// Initial values table may legitimately contain these names
	model_vars.push_back(".RNG.name");
	model_vars.push_back(".RNG.seed");
	model_vars.push_back(".RNG.state");
    }
	
	// Make sure both vectors are sorted to avoid false positive WARNINGs:
	std::sort(model_vars.begin(), model_vars.end());
	std::sort(supplied_vars.begin(), supplied_vars.end());
	
	/*  Test code to check vectors:
	if(data){
		std::cout << "Variables in model:\n";
		std::copy(model_vars.begin(), model_vars.end(),
			  std::ostream_iterator<std::string>(std::cout, ", "));
		std::cout << "\n";

		std::cout << "Supplied vars:\n";
		std::copy(supplied_vars.begin(), supplied_vars.end(),
			  std::ostream_iterator<std::string>(std::cout, ", "));
		std::cout << "\n";
	}  
	*/
	
    std::set_difference(supplied_vars.begin(), supplied_vars.end(),
			model_vars.begin(), model_vars.end(),
			std::inserter(unused_vars, unused_vars.begin()));

    if (!unused_vars.empty()) {
	std::cerr << "\nWARNING: Unused variable(s) in ";
	if (data) {
	    std::cerr << "data table:\n";
	}
	else {
	    std::cerr << "initial value table:\n";
	}
	std::copy(unused_vars.begin(), unused_vars.end(),
		  std::ostream_iterator<std::string>(std::cerr, "\n"));
	std::cerr << "\n";
    }

}

std::string ExpandFileName(char const *s)
{
    if(s[0] != '~') return s;
    std::string name = s;
    if(name.size() > 1 && s[1] != '/') return s;

    char const *p = getenv("HOME");
    if (p) {
	std::string UserHOME = p;
	if (!UserHOME.empty()) {
	    if (name.size() == 1) 
		return UserHOME;
	    else
		return UserHOME + name.substr(1);
	}
    }
    return name;
}


void doSystem(std::string const *command)
{
    std::system(command->c_str());
}

void listModules(){
  std::vector<std::string> mods = jags::Console::listModules();
  std::cout << "Modules:" << std::endl;
  for (unsigned int i = 0; i < mods.size(); ++i) {
      std::cout << "  " << mods[i] << std::endl;
  }
}


void listFactories(jags::FactoryType type)
{
    std::vector<std::pair<std::string, bool> > faclist = 
	jags::Console::listFactories(type);

    std::vector<std::pair<std::string, bool> >::const_iterator p;
    unsigned int max_strlen = 0;
    for (p = faclist.begin(); p != faclist.end(); ++p) {
	if (p->first.length() > max_strlen)
	    max_strlen = p->first.length();
    }
    if (max_strlen < 4)
	max_strlen = 4;

    //Header
    std::cout << "Name";
    for (int i = max_strlen - 4; i >=0; --i) {
	std::cout << " ";
    }
    std::cout << "Status\n";

    //Body
    for (p = faclist.begin(); p != faclist.end(); ++p) {
	std::cout << p->first << " ";
	for (int i = max_strlen - p->first.length(); i >= 0; --i) {
	    std::cout << " ";
	}
	if (p->second) {
	    std::cout << "on";
	}
	else {
	    std::cout << "off";
	}
	std::cout << "\n";
    }
}

void setFactory(std::string const &name, jags::FactoryType type, 
		std::string const &status)
{
    if (status == "on") {
	jags::Console::setFactoryActive(name, type, true);
    }
    else if (status == "off") {
	jags::Console::setFactoryActive(name, type, false);
    }
    else {
	std::cout << "status should be \"on\" or \"off\"";
    }
}

void setSeed(unsigned int seed)
{
    if (seed == 0) {
	std::cout << "seed must be non-zero";
    }
    else {
	jags::Console::setRNGSeed(seed);
    }
}
	    
bool Jtry(bool ok)
{
    if (!ok && !interactive) 
	exit(1);
    else
	return ok;
}

bool Jtry_dump(bool ok)
{
	// Allows doDump to work as described in the manual:
	if (!ok) {
	errordump();
	console->clearModel();
    if (!interactive) 
	exit(1);
	}
	return ok;
}
