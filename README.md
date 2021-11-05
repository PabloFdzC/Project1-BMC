# Project1-BMC
Student:
Jose Pablo Fernández Cubillo - 2019047740

First of all you need to have installed Glade, in
order to compile the project, this can be done with
the command:

sudo apt install glade

Then there is a makefile which compiles the project.
It makes sure that all files and glade is added.

To run the program one goes to the command line and 
writes something like:

./chromosome_maps

Then a window will appear. The user can choose a file
that already exists or proceed to insert the number of
genes needed, this number cannot be greater than 7 nor
lower than 2. After confirming this number the user can
fill the matrix with probabilities, these need to be
numbers lower than 0.5 and greater than 0.0, as soon as
a valid number is inserted the possible chromosomes are
displayed below.

The user can move the chromosomes painted on screen by
dragging the mouse over the drawing section