# tukey-morphisms

TukeyMorphismsBetweenFiniteRelations.R contains code to create (finite) binary relations, represented as either matrices or bipartite graphs, and to determine whether or not a morphism exists between any two such relations. All other files in this repository are outputs from that code. 

The first section of TukeyMorphismsBetweenFiniteRelations.R defines several functions which are useful for generating relations of a given size, calculating the dominating number of a relation, reducing a relation to its skeleton bimorphic form, and finding Tukey morphisms between two relations. 

The second section generates relations up to a certain size. In this case, relations that have up to 6 points in each underlying set. IT IS HIGHLY RECOMMENDED TO NOT RUN ANYTHING GREATER THAN UP TO 5 POINTS ON YOUR LOCAL MACHINE. For context, in the 5x5 case (5 points in each underlying set) there are 324,632 relations. For the 6x6 case, there are 109,453,344 relations. This made my computer very sad. I utilized a High Performance Computing system provided by Boise State University to run this code and produce the associated outputs. (It took approximately 30 hours to run.)

The third section reduces each relation to its skeleton bimorphic form. For the 5x5 case, this results in 32 unique (up to isomorphism) relations. For the 6x6 case, that number is 394.

The fourth section checks for a morphism between each pair of the skeleton relations. 

Finally, the following outputs are produced:

     1. Graphs: contains the graphical representation for the 394 skeleton bimorphic forms.
     2. Skeleton Characteristics_6x6.csv: Lists the dominating number, dual dominating number, and the dual index for each of the 394 skeleton bimorphic forms.
     3. Morphisms_6x6.csv: A 394x394 matrix with a truth value for if there is a morphism from relation i to relation j. 
     4. TukeyMorphism_final.RData: The final R Workspace after running TukeyMorphismsBetweenFiniteRelations.R

Some "fun" facts: 

     - For the 5x5 case, there are 324,632 relations. 
          - This excludes any relation where an element of A- is not related to by something in A+. Since they are trivial with respect to morphism, I don't generate them. 
              - The skeleton bimorphic form for that case, i.e. the empty relation, is added in later.
          - This simplifies to 32 skeleton bimorphic forms.
                - There are 23 unique bimorphic classes (i.e. 9 of the skeleton bimorphic forms are redundant, 7 of them are bimorphic with the 2-ladder.)
          - There are 32x32 = 1,024 pairs to check. The thesis (referenced below) lemmas classify 1,010 (98.63%) of them. 
     - For the 6x6 case, there are 109,453,344 relations. 
          - This simplifies to 394 skeleton bimorphic forms.
               - There are 178 unique bimorphic classes. (i.e. 216 of the skeleton bimorphic forms are redundant, 132 of them are bimorphic with the 2-ladder.)
          - There are 394x394 = 155,236 pairs to check. The thesis lemmas classify 153,066 (98.60%) of them.
     - In the 5x5 case, d(A) = 2 or d(Adual) = 2 for all but the trivial relation (1x1). In the 6x6 case there is an additional relation that does not follow this rule. 
          - #338 is a self-dual relation with a dominating number of 3. 
          - It turns out that there exist self-dual relations for any (finite) dominating number, creating an anti-chain. 

The biggest limitation with the current code is generating all the relations. This takes the majority of the computing time. There might be ways of finding skeleton bimorphic forms more directly, but for now I'm just starting with all possible relations and narrowing it down.

For more background and context, see Barton, Rhett, "Tukey Morphisms Between Finite Binary Relations" (2021). Boise State University Theses and Dissertations. 1868.
https://doi.org/10.18122/td.1868.boisestate
