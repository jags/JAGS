#ifndef GCD_H_
#define GCD_H_

#include <config.h>

#include <atomic>
#include <exception>
#include <mutex>
#include <functional>
#include <cstddef>

#ifdef HAVE_GCD
#include <dispatch/dispatch.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

namespace jags {
    
    // Shared exception state for all worker tasks
    struct ParallelExceptionState {
	std::atomic<bool> has_exception{false};
	std::exception_ptr eptr{nullptr};
	std::mutex m;

 	void capture() noexcept {
	    bool expected = false;
	    if (has_exception.compare_exchange_strong(expected, true,
						      std::memory_order_acq_rel)) {
		std::lock_guard<std::mutex> lock(m);
		eptr = std::current_exception();
	    }
	}
	
	void rethrow_if_any() {
	    if (has_exception.load(std::memory_order_acquire)) {
		std::lock_guard<std::mutex> lock(m);
		std::rethrow_exception(eptr);
	    }
	}
    };
    
#ifdef HAVE_GCD

    /*********************************************************/
    /** GCD Backend for macOS and for Linux with libdispatch */
    /*********************************************************/
    
    inline dispatch_queue_t portable_global_queue() {
#if defined(__APPLE__)
	return dispatch_get_global_queue(QOS_CLASS_USER_INITIATED, 0);
#else
	/* Linux has no QoS API */
	return dispatch_get_global_queue(0, 0);
#endif
    }

    /* Parallel for allowing GCD to manage threads with no
       constraints */
    
    // Context passed to dispatch_apply_f
    struct ApplyContext {
        size_t begin;
        const std::function<void(size_t)> *fn;
	ParallelExceptionState* ex;
    };

    // Trampoline function for dispatch_apply_f
    static void apply_trampoline(void* ctx, size_t i) {
        ApplyContext* c = static_cast<ApplyContext*>(ctx);

	// Skip work if an exception already occurred
        if (c->ex->has_exception.load(std::memory_order_acquire))
            return;
	
        try {
            (*c->fn)(c->begin + i);
        } catch (...) {
            c->ex->capture();
        }
    }

    // Parallel-for using dispatch_apply_f
    inline void parallel_for_apply(size_t begin,
				   size_t end,
				   const std::function<void(size_t)>& fn)
    {
	ParallelExceptionState ex;

	ApplyContext ctx{begin, &fn, &ex};
	size_t count = end - begin;

	dispatch_apply_f(count,
			 portable_global_queue(),
			 &ctx,
			 apply_trampoline);

	// Rethrow on caller thread
	ex.rethrow_if_any();
    }

    /* Parallel for with GCD but limiting concurrency */

    // Context passed to each dispatched task
    struct LimitedContext {
	size_t index;
	const std::function<void(size_t)>* fn;
	dispatch_semaphore_t sem;
	dispatch_group_t group;
	ParallelExceptionState* ex;
    };

    // Trampoline for dispatch_async_f
    static void limited_trampoline(void* p) {
	auto* c = static_cast<LimitedContext*>(p);

	// Skip work if an exception has already occurred
	if (!c->ex->has_exception.load(std::memory_order_acquire)) {
	    try {
		(*c->fn)(c->index);
	    } catch (...) {
		c->ex->capture();
	    }
	}

	dispatch_semaphore_signal(c->sem);
	dispatch_group_leave(c->group);
	delete c;
    }

    // Parallel-for with concurrency limit + exception safety
    inline void parallel_for_limited(size_t begin,
				     size_t end,
				     size_t max_concurrency,
				     const std::function<void(size_t)>& fn)
    {
	dispatch_queue_t q = portable_global_queue();

	if (max_concurrency == 0)
	    max_concurrency = 1;

	dispatch_semaphore_t sem =
	    dispatch_semaphore_create(max_concurrency);

	dispatch_group_t group =
	    dispatch_group_create();

	ParallelExceptionState ex;

	for (size_t i = begin; i < end; ++i) {
	    dispatch_semaphore_wait(sem, DISPATCH_TIME_FOREVER);

	    auto* ctx = new LimitedContext{
		i,
		&fn,
		sem,
		group,
		&ex
	    };

	    dispatch_group_enter(group);
	    dispatch_async_f(q, ctx, limited_trampoline);
	}

	// Wait for all tasks to finish
	dispatch_group_wait(group, DISPATCH_TIME_FOREVER);

	dispatch_release(group);
	dispatch_release(sem);

	// Rethrow the first exception (if any)
	ex.rethrow_if_any();
    }

    /* Convenience wrapper allowing to choose between limited and unlimited
       concurrency. Here nthread = 0 corresponds to no unlimited */
    
    inline void parallel_for(size_t begin,
			     size_t end,
			     const std::function<void(size_t)>& fn,
			     size_t nthread = 0)
    {
	parallel_for_apply(begin, end, fn);
	/*
        if (nthread == 0) {
            parallel_for_apply(begin, end, fn);
        } else {
            parallel_for_limited(begin, end, nthread, fn);
        }
	*/
    }

#else 
    
    /*********************************************************
     OpenMP backend                                       
     If there is no OpenMP then tasks are executed in serial
    *********************************************************/

    inline void parallel_for(size_t begin,
			     size_t end,
			     const std::function<void(size_t)>& fn,
			     size_t nthread = 0)
    {
	ParallelExceptionState ex;
	
	size_t n = end - begin;
	
#ifdef _OPENMP    
#pragma omp parallel for num_threads(nthread > 0 ? nthread : omp_get_max_threads())
#endif
	for (size_t i = 0; i < n; ++i) {

	    // Skip work if another thread already threw
	    if (ex.has_exception.load())
		continue;

	    try {
		fn(begin + i);
	    }
	    catch (...) {
		ex.capture();
	    }
	}

	ex.rethrow_if_any();
    }

#endif /* HAVE_GCD */
    
} // namespace jags

#endif /* GCD_H_ */
