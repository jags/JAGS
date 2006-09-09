#include <config.h>
#include <model/TraceMonitor.h>
#include <graph/Node.h>

#include <stdexcept>
#include <algorithm>

using std::max;
using std::min;
using std::runtime_error;
using std::invalid_argument;
using std::logic_error;

static void align(long &start, long &end, int thin)
{
    /* Ensure that (start % thin == 0) and (end % thin == 0) */

    if (start <= 0)
	start = 1;

    if (end < start)
	end = start;
	
    start = (1 + (start - 1) / thin) * thin;
    end = (end / thin) * thin;
  
    if ((end - start) % thin != 0) {
        throw logic_error("Unable to align start end and thin");
    }
}

TraceMonitor::TraceMonitor(Node const *node, long start, long end, int thin)
{
  _node = node;
  if (thin <= 0) {
    throw invalid_argument("Illegal thinning interval");
  }
  _thin = thin;
  align (start, end, thin);
  _start = start;
  _end = end;
  _size =  (end - start) / thin + 1;
  _values = (double*) calloc(_node->length() * _size, sizeof(double));
  _maxsize = _size;    
  _current = 0;
}

TraceMonitor::TraceMonitor(Node const *node, long start,  int thin)
{
    const long initsize = 128;

    _node = node;
    if (thin <= 0) {
      throw invalid_argument("Illegal thinning interval");
    }
    _thin = thin;
    _size = initsize;
    _start = start;
    _end = start + thin * (_size - 1);
    _values = (double*) calloc (node->length() * _size, sizeof(double));
    _maxsize = -1;
    _current = 0;
}

TraceMonitor::~TraceMonitor ()
{
  free (_values);
}

long TraceMonitor::start() const
{
    return _start;
}

long TraceMonitor::end() const
{
    return _end;
}

long TraceMonitor::thin() const
{
    return _thin;
}

long TraceMonitor::size() const
{
    /* current sample size */
    return _current;
}

long TraceMonitor::size(long start, long end) const
{
  /* sample size between start and end */
  if (end < start)
    return 0;

  if (_current == 0)
    return 0;

  long last_monitored = _start + (_current - 1) * _thin;    
  align (start, end, _thin);

  if (start <= last_monitored && end >= _start) {
    start = max (start, _start);
    end = min (end, last_monitored);
    return (end - start) / _thin + 1;
  }
  else {
    return 0;
  }
}

void TraceMonitor::update(long iteration, unsigned int chain)
{
    if (iteration != _start + _thin * _current) {
	return;
    }

    long node_length = _node->length();

    /* Reallocate vector _values if it is full */
    if (_current == _size && _maxsize == -1) {
      /* Fixme: what do we do if _maxsize != -1 */
      _size = _size * 2;
      _values = (double*) 
	realloc(_values,  node_length * _size * sizeof(double));
      _end = _start + _thin * (_size - 1);
    }
    for (long i = 0; i < node_length; i++) {
      _values[node_length * _current + i] = _node->value(chain)[i];
    }
    _current++;
}

double const *TraceMonitor::values() const
{
  return _values;
}

/*
//FIXME: This belongs somewhere else
void TraceMonitor::dump (ostream &str, long start, long end) const
{
    if (size(start, end) == 0)
	return;

    start = max(start, _start);
    end = min(end, _end);

    if (_node->name().empty()) {
      throw runtime_error("Can't dump nameless node");
    }
    // Write output as an S struct 

    str << "\"" << _node->name() << "\" = structure(c(";
    long node_length = _node->data.length();
    long iter_length = size(start, end);
    for (long j = 0; j < node_length; j++) {
      for (long i = 0; i < iter_length; i++) {
	double v = _values[j + i * node_length];
	if (v == JAGS_NA) {
	  str << "NA";
	}
	else {
	  str << v;
	}
	if (i < iter_length - 1 || j < node_length - 1)
	  str << ",";
      }
    }
    str << "), mcpar = c(" << start << "," << end << "," << _thin << "),";
    if (node_length > 1) {
      // For multivariate nodes, add .Dim attribute 
      str << ".Dim = c(" << iter_length;
      int ndim = _node->data.ndim(false);
      for (long k = 0; k < ndim; k++) {
	str << ",";
	str << _node->data.upper()[k] - _node->data.lower()[k] + 1;
      }
      str << "), ";
    }
    str << "class = \"mcarray\")";
    //str << "Names = \"" << _node->name() << "\")\n";
}
*/

Node const *TraceMonitor::node() const
{
  return _node;
}



