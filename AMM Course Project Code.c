//Computer routine for Structural Analysis.

/* The layout of the program will be as follows -
    -> Ask the user the type of analysis req. Use the Banking program to have different functions for different types of analysis, for instance,
    simple linear, with multi-point constraint, temperature effects, torsion etc. But make the one for this analysis first and only extend it further.

    -> Ask the number of nodes, n. NO. of Elements in n-1. Then ask the coordinates of the nodes, and using this find the length of the elements. Name the
    elements appropriately. Make the element connectivity table. The global stiffness matrix will depend on this. For finding k, go element by element.

    -> Calculate the value of each element stiffness matrix.
    -> Then directly find the value of the Global Stiffness Matrix and Global Force Vector.

    -> Apply the boundary condition. Here, it is req. to know what type of boundary condition is given and using that the global matrices will be modified.
    Applying the boundary condition will render some of the displacements equal to nil. Don't change the other matrices, multiplication with zero will
    automatically make some rows and columns zero.

    -> Now using some method, preferably Gauss-Seidel Iteration, find the unknowns and give them as the answer.
    -> Have a structure for each element. The object type structure will have, element node i, element node j, length of element, E & A of the element.

    -> Have a separate function for getting the Glo+bal Force vector. This can be got directly by asking the user, atleast for the initial basic case. For
    the more complex case the element force vector has to be evaluated.

    -> The nodal displacements are solved for using Gauss-Siedel Iteration, for which 3 for loops are used. The first loop will iterate and the other two
    loops are required to transverse the 2D arrays. Just use the if condition j!=i to sum the product of K and Q, other than the req. one, and then just
    subtract this sum from F.
*/

#include<stdio.h>
#include<math.h>
struct element  //Structure type 1: element.
{
    int i[3]; //Coordinates for Element node i.
    int j[3]; //Coordinates for Element node j.
    int I; //Global node i.
    int J; //Global node j.
    float len;
    float E;
    float A;
    float stiff; //This is the common term multiplied at the start.
    float l;
    float m;
    float n; //The cosines.
};

struct element read (struct element e, int t)  //Function to read the inputs of each element from the user.
{
    printf ("\nGive the global node I and J of element %d.\n", t);
    scanf ("%d %d", &e.I, &e.J);

    printf ("Give the coordinates of element node i of element %d.\n", t);
    scanf ("%d %d %d", &e.i[0],&e.i[1],&e.i[2]); //Here, the ith x,y,z are for ith element.

    printf ("Now, give the coordinates of element node j of element %d.\n", t);
    scanf ("%d %d %d", &e.j[0], &e.j[1], &e.j[2]); //Thus, the x,y,z coordinates of element node j are stored.

    printf("Give the Young's modulus and Area for element %d\n", t);
    scanf("%f %f", &e.E, &e.A);

    //printf ("%d %d %d %d\n", e.i[0], e.i[1], e.i[2], e.J);

    float L=sqrt( pow(e.j[0]-e.i[0],2) + pow(e.j[1]-e.i[1],2) + pow(e.j[2]-e.i[2],2) ); //0 means x, 1 means y and 2 means z.
    e.len = L; //Length of the element.

    e.stiff = e.E*e.A/e.len; //Magnitude of the stiffness.

    e.l = (e.j[0]-e.i[0])/e.len; //cos x.
    e.m = (e.j[1]-e.i[1])/e.len; //cos y.
    e.n = (e.j[2]-e.i[2])/e.len; //cos z.

    printf ("Length of element %d is %f.\n\n", t, e.len);
    return (e);
}

int n,m; //No. of nodes and no. of elements.

int main ()
{
    int t; //The counting variable.
    printf ("Hey, Tell me the no. of nodes and no. of elements.\n");
    scanf ("%d %d", &n, &m); //No. of Nodes and no. of elements.

    struct element a[m];  //Creates m objects of class element.
    for (t=0;t<=m-1;++t)
    {
        a[t] = read(a[t],t); //Function takes the inputs from the user and stores it in each object. This is element t.
        //printf ("Element %d: %f %f %f\n", a[t].len, a[t].E, a[t].A); //Printing the structure.
    }

    //Getting the element connectivity table.
    printf ("\nThe element connectivity table is:\n Element \t Global node i \t Global node j\n");
    for (t=0;t<=m-1;++t)
    {
        printf (" %d \t\t %d \t\t %d \t\t %f\n", t, a[t].I, a[t].J, a[t].stiff);
    }

    struct matrix //structure type 2: 3nx3n matrix.
    {
        float ke[3*n][3*n]; //The element stiffness matrix.
        float f[3*n]; //The element force vector.
    };

    struct matrix b[m]; //m objects of class matrix.

    printf ("\nHey!!!\n");
    //Element stiffness matrix and force vector.
    int r,s;
    printf ("\n%f\n",a[0].stiff);
    for (t=0;t<=m-1;++t)
    {
        printf ("\n%f\n",a[t].stiff);
        // Stiffness Matrix.
        for (r=0;r<3*n;++r)
        {
            //printf ("%f\n\n",b[t].ke[0][0]);
            for (s=0;s<3*n;++s)
            { b[t].ke[r][s]=0;
            printf ("%f\t",b[t].ke[r][s]); } //Initialize all the values to zero.
            printf ("\n");
            b[t].f[r]=0; //similarly for the element force vector.
        }

        //printf ("\n%f %f *********\n",a[t].stiff,pow(a[t].l,2));
        b[t].ke[(3*a[t].I)-2][(3*a[t].I)-2]=a[t].stiff*pow(a[t].l,2); //Giving some positions values.
        b[t].ke[(3*a[t].I)-2][(3*a[t].I)-1]=a[t].stiff*a[t].l*a[t].m;
        b[t].ke[(3*a[t].I)-2][(3*a[t].I)]=a[t].stiff*a[t].l*a[t].n;
        b[t].ke[(3*a[t].I)-2][(3*a[t].J)-2]=a[t].stiff*pow(a[t].l,2)*(-1);
        b[t].ke[(3*a[t].I)-2][(3*a[t].J)-1]=a[t].stiff*a[t].l*a[t].m*(-1);
        b[t].ke[(3*a[t].I)-2][(3*a[t].J)]=a[t].stiff*a[t].l*a[t].n*(-1);

        b[t].ke[(3*a[t].I)-1][(3*a[t].I)-1]=a[t].stiff*pow(a[t].m,2);
        b[t].ke[(3*a[t].I)-1][(3*a[t].I)]=a[t].stiff*a[t].m*a[t].n;
        b[t].ke[(3*a[t].I)-1][(3*a[t].J)-2]=a[t].stiff*a[t].l*a[t].m*(-1);
        b[t].ke[(3*a[t].I)-1][(3*a[t].J)-1]=a[t].stiff*pow(a[t].m,2)*(-1);
        b[t].ke[(3*a[t].I)-1][(3*a[t].J)]=a[t].stiff*a[t].m*a[t].n*(-1);

        b[t].ke[(3*a[t].I)][(3*a[t].I)]=a[t].stiff*pow(a[t].n,2);
        b[t].ke[(3*a[t].I)][(3*a[t].J)-2]=a[t].stiff*a[t].l*a[t].n*(-1);
        b[t].ke[(3*a[t].I)][(3*a[t].J)-1]=a[t].stiff*a[t].m*a[t].n*(-1);
        b[t].ke[(3*a[t].I)][(3*a[t].J)]=a[t].stiff*pow(a[t].n,2)*(-1);

        b[t].ke[(3*a[t].J)-2][(3*a[t].J)-2]=a[t].stiff*pow(a[t].l,2);
        b[t].ke[(3*a[t].J)-2][(3*a[t].J)-1]=a[t].stiff*a[t].l*a[t].m;
        b[t].ke[(3*a[t].J)-2][(3*a[t].J)]=a[t].stiff*a[t].l*a[t].n;

        b[t].ke[(3*a[t].J)-1][(3*a[t].J)-1]=a[t].stiff*pow(a[t].m,2);
        b[t].ke[(3*a[t].J)-1][(3*a[t].J)]=a[t].stiff*a[t].m*a[t].n;

        b[t].ke[(3*a[t].J)][(3*a[t].J)]=a[t].stiff*pow(a[t].n,2);

        printf ("\n\n");
        for (r=0;r<3*n;++r)
        {
            for (s=0;s<3*n;++s)
            { b[t].ke[r][s]=b[t].ke[s][r];
              printf ("%f\t",b[t].ke[r][s]);
            } //using symmetry.
            printf ("\n");
        }
    }

    //Element Force Vector.
    printf ("\nNow, give the force vector components (x,y,z of node i & j resp.)\n");
    float fxi,fyi,fzi,fxj,fyj,fzj;
    for (t=0;t<=m-1;++t)
    {
        printf ("For element %d\n",t);
        scanf ("%f %f %f %f %f %f",&fxi, &fyi, &fzi, &fxj, &fyj, &fzj);

        b[t].f[(3*a[t].I)-2]=fxi;
        b[t].f[(3*a[t].I)-1]=fyi;
        b[t].f[(3*a[t].I)]=fzi;
        b[t].f[(3*a[t].J)-2]=fxj;
        b[t].f[(3*a[t].J)-1]=fyj;
        b[t].f[(3*a[t].J)]=fzj;
    }

    printf("%f\n\nGlobal stiffness matrix is:\n",b[1].ke[1][1]);
    //Global Stiffness Matrix.
    float K[3*n][3*n];
    for (r=0;r<3*n;++r)
    {
        for (s=0;s<3*n;++s)
        { K[r][s]=0; } //Initialize to zero.
    }

    for (r=0;r<3*n;++r)
    {
        for (s=0;s<3*n;++s)
        {
            for (t=0;t<=m-1;++t) //There are m elements.
            { K[r][s]=K[r][s]+b[t].ke[r][s]; } //Sum the i,j element and store in Global Stiffness.
            printf ("%f\t",K[r][s]);
        }
        printf ("\n");
    }

    //Global Force Vector.
    float F[3*n];
    printf ("%f\n\n",F[0]);
    for (r=0;r<3*n;++r)
    { F[r]=0; } //Initialize to zero.

    printf ("\nThe Global Force Vector is - \n");
    for (r=0;r<3*n;++r)
    {
        for (t=0;t<=m-1;++t)
        { F[r]=F[r]+b[t].f[r]; } //Sum the r element of all matrices and save in r of F.
        printf ("%f\n",F[r]);
    }

    //Global Displacement vector.
    float Q[3*n];
    for (r=0;r<3*n;++r)
    { Q[r]=1; } //Initialize to one.


    //Gauss-Siedel Iteration.
    float sum;
    int N=5; //No. of iterations.
    int c;

    printf ("\n");
    for (c=0;c<N;++c)
    {
        printf ("\nIteration number %d: ",c);
        for (r=0;r<3*n;++r)
        {
            sum=0;
            for (s=0;s<3*n;++s)
            {   if (s!=r)
                { sum = sum + Q[s]*K[r][s]; }
            }

            if (K[r][r]!=0)
            { Q[r]=(F[r]-sum)/K[r][r]; }
            else
            { Q[r]=0; }
            printf ("%f\t",Q[r]);
        }
    }

    //Printing the output.
    printf ("\n\nThe nodal displacements are: \n");
    for (r=0;r<3*n;++r)
    {
        printf("Q%d: %f\n",r+1,Q[r]);
    }
}
//**************************************************************//
