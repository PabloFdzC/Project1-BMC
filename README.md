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
genes needed. After confirming this number the user can
fill the matrix with probabilities, these need to be
numbers lower than 0.5 and greater than 0.0, as soon as
a valid number is inserted the possible chromosomes are
displayed below.

Bad news: I cannot figure out some segmentation faults that
sometimes appear, sometimes I get an error then I would
run the program with the debugger and the error does not
appear. So the most I could do is to set a limit to 3
genes, they do display as you insert the probabilities,
but if you erase them it stops working. Infer probabilities
works, file load and save work, thoug with load I sometimes
encounter error that I have not been able to replicate on the
debugger or even running valgrind to see where it might be.

If you run the file in prueba.txt and then infer probabilities
it works.

If you load prueba.txt and then prueba2.txt the drawings will,
but I have been having troubles with setting the correct values
on the text entries.
