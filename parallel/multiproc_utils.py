import mutiprocessing

######################################################################
# Controllers for parallel execution, one per worker.
# Return when a 'None' job (poison pill) is reached.
######################################################################

class Consumer(multiprocessing.Process):
    
    def __init__(self, task_queue, result_queue, block=True, timeout=None):
        multiprocessing.Process.__init__(self)
        self.task_queue = task_queue
        self.result_queue = result_queue
        self.block = block
        self.timeout = timeout

    def run(self):
        proc_name = self.name
        while True:
            next_task = self.task_queue.get(self.block, self.timeout)
            if next_task is None:
                # Poison pill means we should exit
                break
            self.result_queue.put(next_task())
        return
    
class Task(object):
    def __init__(self, func, args):
        self.func = func
        self.args = args
    def __call__(self):
        return self.func(*self.args)
        

