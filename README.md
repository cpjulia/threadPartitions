Usage:

inside main:

variable `numDocs` is the amount of docs in the input set.
RandContext is the struct that handles the random generation of documents. 
Function `generateRandSet()` takes as argument the number of doc ids to generate, and it uses std::uniform_distribution to generate the numbers, but there are, commented out, other non-uniform forms of generating the doc ids. 
`PartitionGenerator`'s constructor takes as argument the budget value and the input set of doc ids (that can be generated with RandContext struct). 
`PartitionGenerator`'s method `partitionDocIdsForThreads()` takes as argument the number of threads and it's the method that will partition the set of doc ids with the algorithm  that finds the gaps and maps the structure. 
`PartitionGenerator`s parittionDocIdsNaive()` uses the naive approach of equally splitting by sectioning the set in halves (which are not the actual halves because of the id offsets).

Restrictions to values:
-budget: any value starting from 0. If budget = 0, the naive approach will be called.
-numDocs: any value starting from 0 that doesn't overflow uint_64
-threads: 1 <= threads <= 16
