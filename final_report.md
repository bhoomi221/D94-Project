# FINAL REPORT -- Reading Course on Scientific & High-Performance Computing



## Questionnaire

1. From the list of topics covered:
```
  [ **(N)** | **(I)** ]  Scientific Computing – Motivation
  [ **(K)** | **(NA)** ]  Introduction/Review of C++
  [ **(K)** | **(NA)** ]  Modular Programming
  [ **(sn)** | **(si)** ]  Libraries
  [ **(sn)** | **(NA)** ]  Make
  [ **(K)** | **(NA)** ]  Version Control
  [ **(K)** | **(NA)** ]  Debugging
  [ **(K)** | **(NA)** ]  Unit Testing
  [ **(sn)** | **(si)** ]  File I/O: self-describing formats, netCDF, HDF5, ...
  [ **(sn)** | **(I)** ]  Numerics
  [ **(N)** | **(I)** ]  Numerical Linear Algebra
  [ **(N)** | **(NA)** ]  ODEs
  [ **(N)** | **(I)** ]  Molecular Dynamics, N-body simulations
  [ **(N)** | **(NA)** ]  PDEs
  [ **(N)** | **(NA)** ]  FFT
  [ **(N)** | **(NA)** ]  Random numbers
  [ **(N)** | **(NA)** ]  Monte Carlo
  [ **(N)** | **(I)** ]  Profiling
  [ **(sn)** | **(I)** ]  Optimization
  [ **(N)** | **(I)** ]  Supercomputing/HPC
  [ **(N)** | **(I)** ]  Shared-memory programming, OMP
  [ **(N)** | **(I)** ]  Distributed-memory programming, MPI
  [ **(N)** |  ]  Other parallel approaches, technologies
  [  |  ]
  [  |  ]
  [  |  ]
  [  |  ]
  [  |  ]
```
**Please add any topic(s) that we may have covered and is not included in the list above.**

1.a) Please describe which topics, from the ones covered in the material were:

   - new **(N)**,
   - not new but some elements were novel **(sn)**,
   - already known **(K)**, i.e. not new.


1.b) Also indicate whether these were:

   - interesting **(I)**,
   - somehow interesting **(si)**,
   - not quite appealing **(NA)**.



2. Please indicate which topics you consider should have been covered with *more* or *less* details?
I think since I was doing the project with MPI it would've been nice to have it covered in a bit more detail with some explanation of other functions that might be useful to use and covering some more examples of how MPI is used.
I also think numerical linear algebra and its uses for AI and machine learning might've been fun to learn about since we are computer science students.


<div style="page-break-after: always; visibility: hidden"> 
\pagebreak 
</div>


10. Please provide any comments about the provided *lectures notes* on
"Practical Scientific Computing" by R.Van Zon and M.Ponce.

<div style="page-break-after: always; visibility: hidden"> 
\pagebreak 
</div>


11. If this would have been a *curricular* course offered to CMS students,
please provide comments, recommendations or suggestions in how to improve,
change, and/or adjust the material covered.

The modular programming, debugging and testing was just reviewed since we have done it in a lot of other classes but the libraries and linking them was good to keep. I think the topics were chosen well. I just  wish there was a bit more depth to them. I felt like 1 class or 2 hours was not enough for parallel programming. Maybe 1 more hour would've been better for both (MPI and OMP). I think in terms of the math topics many people would not be interested in that if it was heavily focussed on the uses of HPC for AI or another more computer science topic like how the supercomputer works it would be more interesting to computer science students. I also liked the nvidia report because it made me realise just how cool HPC can be and what it's already doing for us and it was a very small and easy assignment so there was no pressure. 



11.i Which topics/courses would you say would have been useful to know before
in order to take the best advantage of this course?

B09 - all the knowledge about C and make files
B07- modular programming 

Knowing about debugging is useful as well. 



<div style="page-break-after: always; visibility: hidden"> 
\pagebreak 
</div>


# Report

i. Write a report describing the exercises that you solved,
including which parts you found the most challenging ones for each of them and
how did you tackle them.

Please indicate the location (URL) of the repository where you worked on these
problems.

I didn't solve any exercises I did however do a project on SPH. I found code that did the SPH math in C then I converted it to C++ and then made it compatible for OMP and MPI. I added spatial binning to the code to make the neighbour search easier then I parallelized the code. I was able to add the spatial binning with some issues but I did eventually get there and adding OMP to it was quite easy in comparison. The hardest part was adding MPI. I had a lot of issues with it, for example I had made my serial code send pointers to save memory which I did not know was not allowed for MPI which set me back quite a bit. The hardest part about adding MPI was remembering which process has access to what information because a slight oversight in that aspect would affect all the calculations. So what I had to do was each time there was a bug I would output the information that each process has, to see if there were any differences and if there was I knew I wasn't sharing information that the other processes needed.


<div style="page-break-after: always; visibility: hidden"> 
\pagebreak 
</div>



ii. Write a short report about what you consider is the most interesting
aspects of Scientific Computing and High-Performance Computing.
Also include comments on how you consider (if there is any of) these topics
that could result to be useful for your career path.

These disciplines encompass a wide range of methods, tools, and techniques designed to solve complex computational problems efficiently. I think the simulation and modelling, especially the digital twins are very interesting because the possibilities of what we can do with them is endless. Also the practical uses for these programs can change the world especially when it comes to the climate crisis. Also its implications for quantum computing is very interesting although we don't have quantum computers already, being able to implement some techniques is ground breaking and I'm sure this will eventually lead to a breakthrough in quantum computers. I think the main use for this in my career path will just be scalability, in HPC we are always looking at a large scale and with information data becoming more and more important I think it will help me think about how to deal with such large sets of data or how to keep my work scalable for the future.

<div style="page-break-after: always; visibility: hidden">
\pagebreak
</div>

iii. Consider that you are given a code which is required to be optimized.
The code can take dataset and the size of the data sets affects the different
metrics of the performance (i.e. running time, memory utilization, IO, communication, etc.)
a) What steps would you take in order to study and draft a performance analysis report of this code?

So first you would need to study the code and see what the purpose of the code is on a high level. Then you can decide what metrics you want to focus on for the execution of this code like time vs memory. Then you need to decide what tools you will use to measure these parameters so you can set them up. For example if you want to use an instrumenting method then you will need to go and change the code and decide which points to add the markers. Studying the code from beforehand would definitely help with this. Then you can use this information to identify bottlenecks in the code and provide possible optimizations.

<div style="page-break-after: always; visibility: hidden">
\pagebreak
</div>

iv. You are in the process of submitting an allocation request for compute usage
to one of the national systems, for which you must submit an scaling analysis
of the jobs you will run.
For details on this process, see
	https://alliancecan.ca/en/services/advanced-research-computing/research-portal/resource-allocation-competitions/rac-frequently-asked-questions

a) What type of information should you submit?
b) Should you also submit a strong and weak scaling? Justify.

a) Your request should include:
- A description of your code, the problem it is solving
- The resource requirements along with an explanations of why those resources are important to your work
- Libraries, list any of the external libraries or other software you are using to ensure that its compatible with the system
- Finally submit the result you expect to get from this and the applications of the results


b) Yes, you should submit both of them because the strong scaling will show how efficiently your code uses the additional resources. SO you can see if the speed up is worth the increase in resources. The weak scaling will then show how efficiently it can handle it when the problem size is bigger but there are more resources allocated to it.
These analyses will help the reviewers understand the efficiency and effectiveness of your computational workloads and the benefits of allocating more resources.





