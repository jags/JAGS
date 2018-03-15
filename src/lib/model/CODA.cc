#include <config.h>
#include "CODA.h"

#include <model/Monitor.h>
#include <util/nainf.h>
#include <util/dim.h>

#include <fstream>
#include <sstream>
#include <vector>

using std::list;
using std::vector;
using std::string;
using std::ostringstream;
using std::ofstream;
using std::ostream;

namespace jags {

static void writeDouble(double x, ostream &out)
{
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
    else {
	out << x;
    }
}

    static vector<bool> missingValues(MonitorControl const &control,
				      unsigned int nchain)
    {
	/* 
	   Returns a boolean vector that indicates which variables
	   have at least one missing value in at least one chain.
	   Such variables are omitted when writing the index and
	   output files.
	*/

	Monitor const *monitor = control.monitor();
	unsigned int nvar = product(monitor->dim());
	
	vector<bool> ans(nvar, false);
	for (unsigned int ch = 0; ch < nchain; ++ch) {
	    vector<double> const &y = monitor->value(ch);
	    for (unsigned int v = 0; v < nvar; ++v) {
		if (ans[v]) continue;
		if (monitor->poolIterations()) {
		    if (y[v] == JAGS_NA) {
			ans[v] = true;
		    }
		}
		else {
		    for (unsigned int k = 0; k < control.niter(); ++k) {
			if (y[k * nvar + v] == JAGS_NA) {
			    ans[v] = true;
			    break;
			}
		    }
		}
	    }
	}
	return ans;
    }
    
static void WriteIndex(MonitorControl const &control,
		       vector<bool> const &missing,
		       ofstream &index, unsigned int &lineno)
{
    /* 
       Writes to a coda index file and augments the current line
       number in the corresponding coda output file(s)
    */
    Monitor const *monitor = control.monitor();
    if (monitor->poolIterations()) {
	return;
    }

    unsigned int nvar = product(monitor->dim());
    vector<string> const &enames = monitor->elementNames();
    for (unsigned int v = 0; v < nvar; ++v) {
	if (missing[v]) continue;
	index << enames[v] << " " << lineno + 1 << " "
	      << lineno + control.niter() << '\n';
	lineno += control.niter();
    }
}

//Write output file
static void WriteOutput(MonitorControl const &control, int chain,
			vector<bool> const &missing,
			ofstream &output)
{
    Monitor const *monitor = control.monitor();
    if (monitor->poolIterations()) {
	return;
    }
    
    vector<double> const &y = monitor->value(chain);
    unsigned int nvar = product(monitor->dim());
    for (unsigned int v = 0; v < nvar; ++v) {
	if (missing[v]) continue;
	unsigned int iter = control.start();
	for (unsigned int k = 0; k < control.niter(); ++k) {
	    output << iter << "  ";
	    writeDouble(y[k * nvar + v], output);
	    output << '\n';
	    iter += control.thin();
	}
    }
}

static void WriteTable(MonitorControl const &control, int chain,
		       vector<bool> const &missing,
		       ofstream &index)
{
    Monitor const *monitor = control.monitor();
    if (!monitor->poolIterations()) {
	return;
    }

    vector<double> const &y = monitor->value(chain);
    vector<string> const &enames = monitor->elementNames();
    
    unsigned int nvar = product(monitor->dim());
    for (unsigned int v = 0; v < nvar; ++v) {
	if (missing[v]) continue;
	index << enames[v] << " ";
	writeDouble(y[v], index);
	index << '\n';
    }
}

static bool AnyMonitors(list<MonitorControl> const &mvec,
			bool pooliter, bool poolchains, string const &type)
{
    /* Check for eligible monitors satisfying poolChains and
     * poolIterations */

    list<MonitorControl>::const_iterator p;
    for (p = mvec.begin(); p != mvec.end(); ++p) {
	if (p->monitor()->poolIterations() == pooliter && 
	    p->monitor()->poolChains() == poolchains &&
		( type == "*" || type == p->monitor()->type() )
		) 
	{
	    return true;
	}
    }
    return false;
}

/* CODA output for monitors that do not pool over chains */
unsigned int CODA(list<MonitorControl> const &mvec, string const &stem,
	 unsigned int nchain, string &warn, string const &type)
{
    /* Check for eligible monitors */
    if (!AnyMonitors(mvec, false, false, type))
	return 0;

    /* Open index file */
    string iname = stem + "index.txt";
    ofstream index(iname.c_str());
    if (!index) {
	string msg = string("Failed to open file ") + iname + "\n";
	warn.append(msg);
	return 0;
    }

    /* Open output files */
    vector<ofstream*> output;
    for (unsigned int n = 0; n < nchain; ++n) {
	ostringstream outstream;
	outstream << stem << "chain" << n + 1 << ".txt";
	string oname = outstream.str();
	ofstream *out = new ofstream(oname.c_str());
        if (out) {
            output.push_back(out);
        }
	else {
	    //In case of error, close opened files and return
	    index.close();
	    while(!output.empty()) {
		output.back()->close();
		delete output.back();
		output.pop_back();
	    }
	    string msg = string("Failed to open file ") + oname + "\n";
	    warn.append(msg);
		return 0;
	}
    }
    
    unsigned int lineno = 0;
	unsigned int nwritten = 0;
    list<MonitorControl>::const_iterator p;
    for (p = mvec.begin(); p != mvec.end(); ++p) {
	Monitor const *monitor = p->monitor();
	if (!monitor->poolChains() && !monitor->poolIterations() &&
		( type == "*" || type == monitor->type() ) ) {
	    vector<bool> missing = missingValues(*p, nchain);
	    WriteIndex(*p, missing, index, lineno);
	    for (unsigned int ch = 0; ch < nchain; ++ch) {
		WriteOutput(*p, ch, missing, *output[ch]);
	    }
		nwritten++;
	}
    }

    index.close();
    for (unsigned int i = 0; i < nchain; ++i) {
	output[i]->close();
	delete output[i];
    }
	return nwritten;
}

/* CODA output for monitors that pool over chains */
unsigned int CODA0(list<MonitorControl> const &mvec, string const &stem, string &warn, string const &type)
{
    /* Check for eligible monitors */
    if (!AnyMonitors(mvec, false, true, type))
	return 0;

    /* Open index file */
    string iname = stem + "index0.txt";
    ofstream index(iname.c_str());
    if (!index) {
	string msg = string("Failed to open file ") + iname + "\n";
	warn.append(msg);
	return 0;
    }
    
    /* Open output files */
    string oname = stem + "chain0.txt";
    ofstream output(oname.c_str());
    if (!output) {
	index.close();
	string msg = string("Failed to open file ") + oname + "\n";
	warn.append(msg);
	return 0;
    }
    
    unsigned int lineno = 0;
	unsigned int nwritten = 0;
    list<MonitorControl>::const_iterator p;
    for (p = mvec.begin(); p != mvec.end(); ++p) {
	Monitor const *monitor = p->monitor();
	if (monitor->poolChains() && !monitor->poolIterations() &&
		( type == "*" || type == monitor->type() ) ) {
	    vector<bool> missing = missingValues(*p, 1);
	    WriteIndex(*p, missing, index, lineno);
	    WriteOutput(*p, 0, missing, output);
		nwritten++;
	}
    }
    
    index.close();
    output.close();
	return nwritten;
}

/* TABLE output for monitors that pool over iterations but not over chains
 */
unsigned int TABLE(list<MonitorControl> const &mvec, string const &stem,
	  unsigned int nchain, string &warn, string const &type)
{
    /* Check for eligible monitors */
    if (!AnyMonitors(mvec, true, false, type))
	return 0;

    /* Open output files */
    vector<ofstream*> output;
    for (unsigned int n = 0; n < nchain; ++n) {
	ostringstream outstream;
	outstream << stem << "table" << n + 1 << ".txt";
	string oname = outstream.str();
	ofstream *out = new ofstream(oname.c_str());
        if (out) {
            output.push_back(out);
        }
	else {
	    //In case of error, close opened files and return
	    while(!output.empty()) {
		output.back()->close();
		delete output.back();
		output.pop_back();
	    }
	    string msg = string("Failed to open file ") + oname + "\n";
	    warn.append(msg);
	    return 0;
	}
    }
    
	unsigned int nwritten = 0;
    list<MonitorControl>::const_iterator p;
    for (p = mvec.begin(); p != mvec.end(); ++p) {
	Monitor const *monitor = p->monitor();
	if (!monitor->poolChains() && monitor->poolIterations() &&
		( type == "*" || type == monitor->type() ) ) {
	    vector<bool> missing = missingValues(*p, nchain);
	    for (unsigned int ch = 0; ch < nchain; ++ch) {
		WriteTable(*p, ch, missing, *output[ch]);
	    }
		nwritten++;
	}
    }

    for (unsigned int i = 0; i < nchain; ++i) {
	output[i]->close();
	delete output[i];
    }
	return nwritten;
}

/* TABLE output for monitors that pool over chains and iterations */
unsigned int TABLE0(list<MonitorControl> const &mvec, string const &stem, string &warn, string const &type)
{
    /* Check for eligible monitors */
    if (!AnyMonitors(mvec, true, true, type))
	return 0;

    /* Open output file */
    string iname = stem + "table0.txt";
    ofstream output(iname.c_str());
    if (!output) {
	string msg = string("Failed to open file ") + iname + "\n";
	warn.append(msg);
	return 0;
    }
    
	unsigned int nwritten = 0;
    list<MonitorControl>::const_iterator p;
    for (p = mvec.begin(); p != mvec.end(); ++p) {
	Monitor const *monitor = p->monitor();
	if (monitor->poolChains() && monitor->poolIterations() &&
		( type == "*" || type == monitor->type() ) ) {
	    vector<bool> missing = missingValues(*p, 1);
	    WriteTable(*p, 0, missing, output);
		nwritten++;
	}
    }
    
    output.close();
	return nwritten;
}

} //namespace jags
