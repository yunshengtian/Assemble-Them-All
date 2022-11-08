import multiprocessing as mp
from tqdm import tqdm


def parallel_worker(worker, args, queue, proc_idx):
    result = worker(*args)
    queue.put([result, proc_idx])


def parallel_execute(worker, args, num_proc, show_progress=True, desc=None, terminate_func=None):
    '''
    Tool for parallel execution
    '''
    if show_progress:
        pbar = tqdm(total=len(args), desc=desc)

    queue = mp.Queue()
    procs = {}
    n_active_proc = 0

    try:

        # loop over arguments for all processes
        for proc_idx, arg in enumerate(args):

            if num_proc > 1:
                proc = mp.Process(target=parallel_worker, args=(worker, arg, queue, proc_idx))
                proc.start()
                procs[proc_idx] = proc
                n_active_proc += 1

                if n_active_proc >= num_proc: # launch a new process after an existing one finishes
                    result, proc_idx = queue.get()
                    procs.pop(proc_idx)
                    yield result

                    if terminate_func and terminate_func(result): # terminate condition meets
                        for p in procs.values(): # terminate all running processes
                            p.terminate()
                        if show_progress:
                            pbar.update(pbar.total - pbar.last_print_n)
                            pbar.close()
                        return
                    
                    n_active_proc -= 1

                    if show_progress:
                        pbar.update(1)
            else:
                parallel_worker(worker, arg, queue, proc_idx) # no need to use mp.Process when serial
                result, _ = queue.get()
                yield result

                if terminate_func and terminate_func(result): # terminate condition meets
                    if show_progress:
                        pbar.update(pbar.total - pbar.last_print_n)
                        pbar.close()
                    return

                if show_progress:
                    pbar.update(1)

        for _ in range(n_active_proc): # wait for existing processes to finish
            result, proc_idx = queue.get()
            procs.pop(proc_idx)
            yield result

            if terminate_func and terminate_func(result): # terminate condition meets
                for p in procs.values(): # terminate all running processes
                    p.terminate()
                if show_progress:
                    pbar.update(pbar.total - pbar.last_print_n)
                    pbar.close()
                return

            if show_progress:
                pbar.update(1)

    except:
        for proc in procs.values():
            proc.terminate()

    if show_progress:
        pbar.close()
