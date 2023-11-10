from multiprocessing import Process, JoinableQueue
import time

def fun1(a):
    print(f"pass fun1:{a}")
    time.sleep(2)


def fun2(a):
    print(f"pass fun2:{a}")
    time.sleep(2)


def fun3(a):
    print(f"pass fun3:{a}")
    time.sleep(2)


def add(queue):
    for f in ['fun1', 'fun2', 'fun3']:
        task = (f, 1)
        print(f"add task {f}")
        queue.put(task)
        # time.sleep(5)
    queue.put(None)

def consume(queue):
    while True:
        task = queue.get()
        if task is None:
            break
        p = Process(
            target=globals()[task[0]],
            args=(task[1:],)
        )
        p.start()
        queue.task_done()
    queue.task_done()



if __name__ == "__main__":
    queue = JoinableQueue()
    p1 = Process(target=add, args=(queue,))
    p1.start()

    p2 = Process(target=consume, args=(queue,))
    p2.start()
    p2.join()
    # queue.join()
    print('Done')

