# Task-Based Programming with TBB

## Benefits of Task-Based Programming

- Matches parallelism to available resources
- Provides faster task startup and shutdown
- Allows more efficient evaluation order
- Improves load balancing
- Enables higher–level thinking

## Detailed Benefits

- Threads map onto the physical threads of the hardware for high efficiency
- Tasks in TBB are lighter weight than logical threads
- Tasks in TBB are efficient because the scheduler is unfair, it sacrifices fairness for efficiency
- The scheduler does load balancing, distributing work evenly across threads
- Task-based thinking allows focus on logical dependences between tasks, letting the scheduler handle efficiency

**Tip:**
> Design programs to create many more tasks than there are threads. Let the task scheduler map tasks to threads.

## When Task-Based Programming Is Inappropriate

- Not suitable for tasks that block frequently (e.g., waiting for I/O or mutexes for long periods)
- For blocking tasks, use full-blown threads

## How Task Scheduler Works

- Enables many threads by creating enough jobs
- Preserves data locality to make single thread execution more efficient
- Minimizes both memory demands and cross-thread communication

**Strategies:**

- Balance between depth-first and breadth-first execution
- Depth-first is better for sequential execution due to cache efficiency and minimizing space
- Each thread has a deque of tasks ready to run
- Tasks are executed based on a set of rules that ensure efficient utilization

## Task Scheduler Bypass

- Directly specify the next task to run to avoid unnecessary deque operations
- Guarantees that the task is executed by the current thread and not by any other thread

**Usage:**
> Available only with the preview feature of `onepai::tbb::task_group`

# Design Patterns in oneAPI Threading Building Blocks (oneTBB)

Each pattern is described in the following format:

- Problem
- Context
- Forces
- Solution
- Example

The patterns listed are:

- Agglomeration
- Elementwise
- Odd-Even Communication
- Wavefront
- Reduction
- Divide and Conquer
- GUI Thread
- Non-Preemptive Priorities
- Lazy Initialization
- Local Serializer
- Fenced Data Transfer
- Reference Counting

## Agglomeration

### Problem

Parallelism at a very fine grain level leads to synchronization overhead.

### Context

Algorithms with very fine-grain parallelism can be overwhelmed by the synchronization overhead between threads.

### Forces

- Computations are small and can be done in parallel.
- The parallelism is for performance, not semantics.

### Solution

Group computations into blocks and evaluate computations within a block serially. The block size should be large enough
to amortize parallel overhead, but not so large that it limits parallelism or load balancing.

### Example

oneTBB loop templates like `oneapi::tbb::parallel_for` support automatic agglomeration. Cache effects and boundary to
interior ratio effects should be considered. Also, blocks containing long contiguous subsets of data may better enable
vectorization.

For recursive computations, the solution is to treat subtrees as groups. An example might be a recursive sort, which
solves sub-problems in parallel only if they are above a certain threshold size.

### Reference

[Designing and Building Parallel Programs](http://www.mcs.anl.gov/~itf/dbpp) by Ian Foster introduced the term "
agglomeration". The book discusses a four-step PCAM design method:

- Partitioning: break the program into the smallest tasks possible.
- Communication: figure out required communication between tasks. In oneTBB, communication is usually cache line
  transfers.
- Agglomeration: combine tasks into larger tasks.
- Mapping: map tasks onto processors. The oneTBB task scheduler handles this step.