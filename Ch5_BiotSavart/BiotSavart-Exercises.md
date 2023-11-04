# Calculating the magnetic field with the Biot-Savart Law

Developed by J. D. McDonnell

In this set of exercises, the student will implement code for the Biot-Savart law to compute the magnetic field at any point in space due to a square-shaped loop.  

Upon calculating the magnetic field at many points in space, the student will prepare a vector plot of the magnetic field.

## Exercises
### Exercise 1: Calculate the magnetic field due to a straight segment of wire

Consider a straight segment of wire $1$ unit long.  Place one end at $(0.5, 0.5, 0.0)$, and the other end at $(-0.5, 0.5, 0.0)$.  Let a steady current of $1$A flow through this wire, from the first end towards the second end.  (Such a current cannot exist physically, but this is still a good first-step towards the goal of these Exercises.)

1. Use the Biot-Savart law to analytically calculate the magnetic field at the origin.  This will serve as a check for the numerical method.  
2. Describe in words (or pseudocode) a procedure to *numerically* calculate the magnetic field at the origin with the Biot-Savart law.  
3. Implement the code to numerically calculate the magnetic field at the origin with the Biot-Savart law.  Make sure that your numerical answer agrees with the analytical answer.  

### Exercise 2: Calculate the magnetic field at the center of a square loop

Consider a square loop of wire, sitting in the $xy$-plane. There is a current of $1$A flowing through the wire. The goal of these exercises is to calculate the magnetic field that results from this current configuration.

1. Use the Biot-Savart law to analytically calculate the magnetic field at the center of the square loop - assume the loop has sides of length $1$ unit for simplicity. This will serve as a check for the numerical method.  **Note**: For both this analytical calculation and your numerical calculation below, think of this square loop of wire as four segments of straight wire connected to each other.  You can then build from your work in Exercise 1.  
2. Describe in words (or pseudocode) a procedure to *numerically* calculate the magnetic field at the center of the loop with the Biot-Savart law.  
3. Implement the code to numerically calculate the magnetic field at the center of the loop with the Biot-Savart law.  Make sure that your numerical answer agrees with the analytical answer.

### Exercise 3: Calculate the magnetic field of a square loop at any point in space

Now that you have validated your numerical approach for the magnetic field at the center of the loop, your new task is to *generalize* your approach in order to calculate the magnetic field at *any* point in space $\vec{r}$.  In particular, you will calculate the magnetic field at many grid points in the $yz$-plane, and from there you will be able to visualize the magnetic field with a vector plot.

1. Describe in words (or pseudocode) the modifications you will need to make to your previous procedure.  Discuss the way that you will calculate the magnetic field at many points in the $yz$-plane.
2. Implement the code to calculate the magnetic field at many points in the $yz$-plane.
3. From the magnetic field that you have calculated, produce a vector plot. Describe the key features that you see. Does the image match your expectations?  Is there anything that surprises you? Where is the magnetic field the strongest, and where is it weakest? How can you tell from the field lines?


### Extension: Calculate the magnetic field for any wire shape

The first two exercises focus on calculating the magnetic field due to a square loop of current-carrying wire. It is natural to extend this exercise to calculate the magnetic field due to *other* interesting shapes of curent-carrying wire. Some suggested shapes:

  - A circular loop of wire. In this case, you will be able to perform an analytical verification of the magnetic field at the center of the loop again.  
  - A helix of wire. If you make the helix long enough and with tight "windings", you can approximate a finite-length solenoid.


