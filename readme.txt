
This code simulates the evolution of influenze... 


The central entity is the class CStrain which modulates all 
the info related to a strain, such as number of people infected 
by that (N), fitness and a vector which contains the pointer to 
the neighbouring nodes in the tree including its father as its 
first element (except when it has no father). See strain.h for 
details of the interface. 

The list strains contains all the alive strains. Once a 
strain extinguished it is removed from the list and it is signaled 
as dead and it (deallocatable) memory is freed to save some memory, but 
the back bone is kept specially because we want to preserve the tree.

Things to be checked and improved:
- Virtually every thing need to be checked
- The solver needs to be dramatically improved
- Check the results against an analytically known case
- Look for places where the code can be optimized
