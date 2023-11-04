# Charges in a conductor and Gauss's Law

Developed by Larry Engelhardt

In this activity, students simulate and visualize the motion of excess charges in a conductor (both a sphere and a cube) in three-dimensional space, without any external electric field.  Students also compute the electric field from the conductor and discuss the results in the following three regimes: (1) Far from the conductor, (2) Close to the conductor, and (3) Inside the conductor. 

Note: Student-facing documents that have been used with students in a lab setting are provided under the "Code" tab as "Additional Resources".

## Exercises
**Exercise 1**:  Number of electrons

Throughout these exercises, we will simulate the motion of excess charges in a conductor for a total excess charge of $Q = -5$ $\mu$C.  Given this value of $Q$, calculate the number of excess *electrons* in the conductor.

**Exercise 2**:  Observing the motion of $N = 200$ charges in a conducting sphere

For Exercise 1, you should have calculated a very large number of electrons.  It would be nice to visualize the position of *every* electron in the conductor, but a computer can only do about a billion mathematical operations per second, and if a computer were made much faster, it would melt!  Hence, it is impossible to visualize every electron, even with today’s fastest computers.  Instead, we will treat this net charge as being composed of $N = 200$ particles, each with charge $q = Q/N$, so each of the $N$ charges will actually represent many electrons.

You will be provided with a working computer program that simulates the motion of $N$ charges inside of a conducting sphere of radius $R = 0.1$ m. This simulation could be implemented in any programming language, but an example program is provided here as an online "trinket" implemented using Glowscript/Vpython:

[https://trinket.io/library/trinkets/88da27559b](https://trinket.io/library/trinkets/88da27559b)

*TRINKET TIP: In order to clearly see the contents of all of the windows in this trinket, you should display the trinket in Fullscreen mode by clicking the menu in the top-left corner of the trinket and selecting “Fullscreen”. Then click the "Run" button at the top to execute the simulation. You can control-drag to look from different angles, and alt-drag to zoom in and out.*

In this program, the charges start out being distributed throughout the volume of the sphere with random initial positions.  Then the charges exert forces on one another and move.  What happens to the charges?  Where do the charges end up? Why does this make sense?

After the charges have reached equilibrium, take a screenshot of the arrangement of the charges, and briefly describe the arrangement of charges.


**Exercise 3**:  Writing a function to compute the electric field

Inside of the program that you just executed, you are provided with an **incomplete** function titled `computeEField`. Complete this function so that the function will compute the net electric field at point P in units of N/$\mu$C. (Notes: The reason for converting to units of N/$\mu$C is that some of the numbers that you will calculate would be quite large if you left them in units of N/C, so converting to N/$\mu$C will make the numbers easier to read. The detailed syntax of this function will depend on which programming language you use, but the following tip will apply regardless of implementation language.)

**Tip**:  Look at the function `computeForces` that is defined just above `computeEField` and copy and paste code from `computeForces` into `computeEField`.  (Computing the net electric field is very similar to computing the net forces, but computing the net electric field is ***simpler*** since you don’t need to go through every *pair* of particles.)

Remember, whenever you write a computer program, it won’t work right away.  That is okay!  Try things, be patient, and ask for help as needed.

Once you have completed your code, also briefly describe what your code does.

**Exercise 4**: Validating the E-field far from a conducting sphere: $r = 0.5$ m, $R = 0.1$ m

Using Coulomb’s Law, calculate the electric field that you would expect to observe at a distance of $r = 0.5$ meters away from a charge $Q = -5$ $\mu$C.  

Execute your program using $N = 200$ charges, and use your program to compute the electric field at a point P located along the x axis, a distance $r = 0.5$ meters away from the center of the conducting sphere: $\vec{P}$ = {0.5 m, 0, 0}

Verify that the electric field that you compute from the distribution of $N = 200$ charges is indeed consistent with Coulomb’s Law for $r = 0.5$ meters.  If the answers are not consistent, then something is wrong that you will need to fix.

**Exercise 5**:  E-field close to a conducting sphere:  $r = 0.15$ m and $R = 0.10$ m

In your program, change the position of point P to 
$\vec{P}$ = {0.15 m, 0, 0}
which is relatively close to the spherical conductor of radius $R = 0.1$ m.

Execute your program using this value of P, and compare your results with the value predicted by Gauss’ Law for $r = 0.15$ m.  Why/how is Gauss’s Law relevant here?

**Exercise 6**: E-field *inside* of a conducting sphere:  $r = 0.05$ m and $R = 0.10$ m

Now consider a point that is *inside* the conducting sphere (for example, $r = 0.05$ m and $R = 0.10$ m). How do you think the electric field at this point will compare the electric field that you computed in Exercise 5?

In your program, change the position of point P to $\vec{P}$ = {0.05 m, 0, 0} which is inside the spherical conductor of radius $R = 0.1$ m.

Execute your program using this value of P.  What do you observe for the electric field?  In particular, how does the value of $E$ change while the $N = 200$ charges approach equilibrium?  What value do you observe for $E$? Is this what you expected?  How is this related to Gauss’ Law? Feel free to experiment with this, repeating the computation at different points, P, inside of the conducting sphere.

**Exercise 7**:  Observing the motion of $N = 8$ charges in a conducting cube

You will again simulate the motion of $N$ charges in a conductor, but now the conductor will be shaped as a ***cube*** instead of being shaped as a ***sphere***. You will again be provided with a working computer program which could be implemented in any programming language, but an example program is provided here as an online trinket implemented using Glowscript/Vpython:

[https://trinket.io/library/trinkets/e8720c9884](https://trinket.io/library/trinkets/e8720c9884)

Execute this program using $N = 8$ charges, and observe what the charges do.  Repeat the simulation a few times.  What happens, and why does is make sense?  

After the charges have reached equilibrium, take a screenshot of the arrangement of the charges, and briefly describe the arrangement of charges.

**Exercise 8**:  Observing the motion of $N = 20$ charges in a conducting cube

Repeat Exercise 7, this time using $N = 20$ charges inside the conducting cube.  Again, take a screenshot of the arrangement of the charges, and briefly describe the arrangement of charges.

**Exercise 9**:  Observing the motion of $N = 200$ charges in a conducting cube

Repeat Exercise 8, this time using $N = 200$ charges inside the conducting cube, and look very closely at the arrangement of the charges.  Is there anywhere that the charges are closer together?  Is there anywhere that the charges are farther apart?  

Again, take a screenshot of the arrangement of the charges, and briefly describe the arrangement of charges.


**Exercise 10**:  Validating the E-field far from a conducting cube: $r = 0.5$ m, $L = 0.1$ m

In Exercise 3 you wrote some code (a "function") to compute the net electric field produced by N charges.  Copy and paste this function into the simulation of the cube, and add a line to call this function. 

Execute your program using $N = 200$ charges with a net charge of $Q = -5$ $\mu$C, and compute the electric field at a point P located along the x axis, a distance $r = 0.5$ meters away from the center of the conducting cube: $\vec{P}$ = {0.5 m, 0, 0}

Compare your results with Coulomb’s Law for a point charge (from **Exercise 4**), and verify that the electric field that you compute from the distribution of $N = 200$ charges is indeed consistent with Coulomb’s Law for $r = 0.5$ meters.  If the answers are not consistent, then something is wrong that you will need to fix.

**Exercise 11**:  E-field close to a conducting cube:  $x = 0.1$ m and $L = 0.1$ m

In your program, change the position of point P to   $\vec{P}$ = {0.1 m, 0, 0}
which is relatively close to the conducting cube.  (The cube is centered at the origin with a length of $L = 0.1$ m, extending from $x = -L/2$ to $x = +L/2$.)

Execute your program using this value of P, and record the value of the electric field produced by the charged cube.  Can you use Gauss’ Law to calculate an analytical value of the electric field at this point P?  Why or why not?

**Exercise 12**:  E-field inside a conducting cube:  $x = 0.01$ m and $L = 0.1$ m

In your program, change the position of point P to   $\vec{P}$ = {0.01 m, 0, 0}
which is 1 cm away from the center of the conducting cube.

Execute your program using this value of P.  What do you observe for the electric field?  In particular, how does the value of E change while the $N = 200$ charges approach equilibrium?  What value do you observe for $E$?  (To get accurate results, you might need to increase $N$ to $N = 500$ or $N = 1000$.)  How are your results related to Gauss’ Law?


