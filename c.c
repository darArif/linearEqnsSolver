/*
   C Program for Solving a System of Linear Equations Through Elementary Row Transformations

   Abbreviations Used:
   REF = Row Echleon Form (Obtained by Applying Gaussian Elimination on a Matrix)
   RREF = Reduced Row Echleon Form  (Obtained by Applying Jordan Elimination on REF)
*/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>

//The Following enum data type will represent if the Entered System of Linear Equations is Homogeneous or Non-Homogeneous
typedef enum 
{
    homogeneous,
    non_homogeneous
} type;

unsigned int lSearch(char *, char); //This Function Returns the Frequency of the Character in the string

void aug_mtrx(float **, float*, int, int); //This Fucntion Converts the Coefficient Matrix into Augmented Matrix

void g_elim(float **, int, int); //This Function Carries Out Gaussian Elimination on the Augmented Matrix

void swap(float *, float *); //This is a Simple Swap Function that Uses Temporary Variable 

int rank_mtrx(float **, int, int); //This Function Returns the Rank of Coefficient Matrix Part of the Augmented Matrix 

void j_elim(float **, int, int); /*This Function Carries Out Jordan Elimination on RREF Created by g_elim() function. 
Also the code for this Jordan Elimination Fucntion has been designed only for an REF with rank=m=n 
where m is the Number of Rows in the Augmented and the Coefficient Matrix and n is the Number of Columns in the Coefficient
Matrix only(the Augmented Matrix has n+1 Columns) */

int main(void)
{
    while (1)
    {
        int m; //Stores the number of Linear Equations
        int n; //Stores the number of Variables
        int ch; //For Removing Unwanted Input Buffer

        int i, j, k;

        printf("\nEnter the no. of Linear Equations : ");
        scanf("%d", &m);
        while ( ( ch = getchar() ) != '\n'); //Removes \n character from Input Buffer

        printf("\nEnter the no. of Variables : ");
        scanf("%d", &n);
        while ( ( ch = getchar() ) != '\n'); //Removes \n character from Input Buffer

        //The following Dynamically Created Matrix M will First Represent the Coefficient Matrix and then the Augmented Matrix of the Entered System  
        float **M = (float **)malloc(m*sizeof(float *));
        for (i=0; i<m; i++)
            M[i] = (float *)malloc((n+1)*sizeof(float));

        //The following Dynamically Created 1-D Array will Store the Names of the Variables
        char *X = (char *)malloc(n*sizeof(char));

        //The following Dynamically Created 1-D Array will Store the Constant Terms of the Entered Equations
        float *B = (float *)malloc(m*sizeof(float));

        //Instructions for Entering the Equations
        printf("\nRead the instructions for entering the equations :\n\n"
           "The standard forms of the Eqns have to be entered one by one in the format : "
           "ax +by +cz =d\nMind the whitespaces!\n"
           "Variables must be single charactered.\n"
           "0 &/ 1 coefficient(s), if any, are necessary.\n"
           "\nFor example, the equations :\n\n"
           "x - 2/5z = 6\n\n2y + 15z = -7\n\n-x + 2z + 3y -8 = 0\n\n"
           "must respectivelly be entered as :\n\n"
           "5x +0y -2z =30\n\n0x +2y +15z =-7\n\n-1x +3y +2z =8\n");

        printf("\nIf you have read the instructions, press Enter twice to continue : ");
        ch = getchar();
        while ( ( ch = getchar() ) != '\n'); //Removes \n character from Input Buffer

        for (i=0; i<m; i++)
        {
            //This Variable Sized Array is for Storing Entered Equations Temporarily
            char strEqn[30*n]; //Since Variable Sized Stack Arrays are Supported in C(from C99)

            char s[30]; //This Array is for Temporarily Storing Tokens Extracted from strEqn String

            lbl6 :
            {
                printf("\nEnter the Eqn %d : ", i+1);
                fgets(strEqn, 30*n, stdin);
            }

            if (lSearch(strEqn, ' ') != n) //Valid Entered Equation has n Single-White-Spaces
                {
                    printf("\nEqn Entered in Invalid Format!\nEnter it again!\n");
                    goto lbl6;
                }

            j = 0;

            char *tok = strtok(strEqn, " "); //Extraction of Tokens Containing the Coefficients and the Variables

            while (tok != NULL)  //For Loading the Extracted Values in Marix M and Array B
            {
                if ( strchr(tok, '=') == NULL )
                {
                    strcpy(s, tok);

                    int len = strlen(s);
                
                    for (k=len; k>(len-2); k--)
                        s[k+1] = s[k];

                    s[len-1] = ' ';

                    if (i == 0)
                        sscanf(s, "%f %c", &M[i][j], &X[j]);

                    else
                        sscanf(s, "%f %*c", &M[i][j]);
                }

                else 
                    sscanf(tok, "%*c %f", &B[i]);

                j++;
                tok = strtok(NULL, " ");
            }
        }
        
        type system = homogeneous;

        for (i=0; i<m; i++)
        {
            if (B[i] != 0) //If Atleast One Constant Term is Non-Zero, the System is Non-Homogeneous
            {
                system = non_homogeneous;
                break;
            }
        }

        int response1; //Stores Response for Detailed or Short Solution

        lbl1 :
        {
            printf("\nFor detailed solution, enter 1\nelse enter 0 : ");
            scanf("%d", &response1);
            while ( ( ch = getchar() ) != '\n');
        }

        if (response1 != 1 && response1 != 0)
        {
            printf("\nInvalid Integer Entered!\nEnter Again!");
            goto lbl1;
        }

        if (response1 == 1)
        {
            printf("\nCoefficient Matrix of the System is :\n{\n");
            for (i=0; i<m; i++) 
            {
                printf(" {    ");
                for (j=0; j<n; j++) printf("%.2f    ", M[i][j]);
                printf("}\n");
            }
            printf("}\n");
        }
        

        aug_mtrx(M, B, m, n);

        if (response1 == 1)
        {
            printf("\nAugmented Matrix of the system is :\n{\n");
            for (i=0; i<m; i++) 
            {
                printf(" {    ");
                for (j=0; j<(n+1); j++) printf("%.2f    ", M[i][j]);
                printf("}\n");
            }
            printf("}\n");
        }
        
        //ref Matrix Stores the REF of the System Formed by Guassian Elimination
        float **ref = (float **)malloc(m*sizeof(float *));
        for (i=0; i<m; i++)
        {
            ref[i] = (float *)malloc((n+1)*sizeof(float));
            for (j=0; j<(n+1); j++)
                ref[i][j] = M[i][j];
        }

        g_elim(ref, m, n+1);

        if (response1 == 1)
        {
            printf("\nRow Echelon Form of the Augmented Matrix is :\n{\n");
            for (i=0; i<m; i++) 
            {
                printf(" {    ");
                for (j=0; j<(n+1); j++) printf("%.2f    ", ref[i][j]);
                printf("}\n");
            }
            printf("}\n");
        }
        
        int rank = rank_mtrx(ref, m, n+1);

        if (response1 == 1)
            printf("\nRank of the Coefficient Matrix is %d\n", rank);

        if (rank == m && rank == n)
        {
            if (system == homogeneous)
            {
                printf("\nThe given homogeneous system of equations has a unique solution, the trivial soluton, i.e., \n\n");
                for (j=0; j<n; j++)
                    printf("%c = 0\n", X[j]);
                printf("\n");

                goto last_lbl;
            }

            else if (system == non_homogeneous)
                goto lbl2;
        }

        else if (rank == n && rank < m)
        {
            if (system == homogeneous)
            {
                printf("\nThe given homogeneous system of equations has a unique solution, the trivial soluton, i.e.,\n\n");
                for (j=0; j<n; j++)
                    printf("%c = 0\n", X[j]);
                printf("\n");

                goto last_lbl;
            }

            else if (system == non_homogeneous)
            {
                goto lbl2;
            }
        }

        else if (rank == m && rank < n)
        {
            if (system == homogeneous)
            {
                printf("\nThe given homogeneous system of equations is consistent with infinite solutions,\n"
                       "both trivial one and non-trivial ones.");
                
                goto last_lbl;
            }

            else if (system == non_homogeneous)
            {
                printf("\nThe given non-homogeneous system of equations is consistent with infinite solutions.");

                goto last_lbl;
            }
        }

        else if (rank < m && rank < n)
        {
            if (system == homogeneous)
            {
                printf("\nThe given homogeneous system of equations is consistent with infinite solutions,\n"
                       "both trivial one and non-trivial ones.");
                
                goto last_lbl;
            }

            else if (system == non_homogeneous)
                goto lbl3;
        }

        lbl2 :
        {
            //rref Matrix Stores the RREF of the System Formed by Jordan Elimination of REF of the System 
            float **rref = (float **)malloc(m*sizeof(float *));

            for (i=0; i<m; i++)
            {
                rref[i] = (float *)malloc((n+1)*sizeof(float));
                for (j=0; j<(n+1); j++)
                    rref[i][j] = ref[i][j];
            }

            j_elim(rref, m, n+1);

            if (response1 == 1)
            {
                printf("\nReduced Row Echelon Form of the Augumented Matrix is :\n{\n");
                for (i=0; i<m; i++) 
                {
                    printf(" {    ");
                    for (j=0; j<(n+1); j++) printf("%.2f    ", rref[i][j]);
                    printf("}\n");
                }
                printf("}\n");
            }
            
            printf("\n\nThe given non-homogeneous system of equations is consistent with a unique solution given as :\n\n");
            for (j=0; j<n; j++)
                printf("%c = %.3f\n", X[j], rref[j][n]);
            printf("Note:- In case of answers in fractions, the answers have been rounded off upto 3 decimal places.\n");


            for (i=0; i<m; i++)
                    free(rref[i]);
                free(rref);

            goto last_lbl;
        }

        lbl3 :
        {
            for (i=0; i<m; i++)
            {
                for (j=0; j<n; j++)
                {
                    if (ref[i][j] != 0)
                        break;
                }

                if (j==n && ref[i][j-1] == 0)
                {
                    if (ref[i][n] != 0)
                    {
                        printf("\nThe given non-homogeneous system of equations is inconsistent and has no solution.\n");
                        goto last_lbl;
                    }
                }
            }

            printf("\nThe given non-homogeneous system of equations is consistent with infininte solutions.\n");
            goto last_lbl;
        }

        last_lbl :
        {
            int response2;

            printf("\n\nFor Solving one more system, enter 1\n"
                   "For exiting, enter 0 :\n");
            scanf("%d", &response2);
            while ( ( ch = getchar() ) != '\n');

            if (response2 == 1)
            { 
                for (i=0; i<m; i++)
                    free(M[i]);
                free(M);

                for (i=0; i<m; i++)
                    free(ref[i]);
                free(ref);

                free(X);

                free(B);

                continue;
            }

            else if (response2 == 0)
            {
                for (i=0; i<m; i++)
                    free(M[i]);
                free(M);

                for (i=0; i<m; i++)
                    free(ref[i]);
                free(ref);

                free(X);

                free(B);

                return 0;
            }

            else
            {
                printf("\nInvalid Integer Entered!\nEnter again!");
                goto last_lbl;
            }
        }
    }
    return 0;
}

void aug_mtrx(float **M, float *B, int m, int n)
{
    int i;

    for (i=0; i<m; i++)
        M[i][n] = B[i];
}

void g_elim(float **ref, int rows, int colms)
{
    int flag;

    float mul, div;

    int i, j, k, l, m, n, o, p, q;

    for (i=0; i<rows; i++)
    {
        j = i;

        lbl1 :
        {
            if (ref[i][j] != 0 ) flag = 1;

            else
            {
                flag = 0;

                k = i + 1;

                while (k < rows)
                {
                    if (ref[k][j] != 0)
                    {
                        for (l=j; l<colms; l++) swap(&ref[i][l], &ref[k][l]);
                           
                        flag = 1;

                        break;
                    }

                    else k++;
                }
            }

            if (flag == 0 && j < colms)
            {
                j++;
                goto lbl1;
            }

            else if (flag == 0 && j >= colms) return;

            else if (flag == 1)
            {
                div = ref[i][j] ;
                for (m=j; m<colms; m++) ref[i][m] = ref[i][m] / div ;
                    
                for (n=i+1; n<rows; n++)
                {   
                    mul = ref[n][j] ;

                    for (o=j; o<colms; o++)   
                        ref[n][o] = ref[n][o] - (ref[i][o])*(mul) ;
                }
            }
        }
    }
    
}
    
void swap(float *p1, float *p2)
{
   float temp = *p1;
   *p1 = *p2;
   *p2 = temp;
}

int rank_mtrx(float **ref, int rows, int colms)
{
    int rank = 0;

    int i, j;

    for (i=0; i<rows; i++)
    {
        for (j=0; j<(colms-1); j++)
        {
            if (ref[i][j] != 0)
            {
                rank++;
                break;
            }
        }
    }

    return rank;
}

void j_elim(float **rref, int rows, int colms)
{
    int j=1, k, l, m, n=1;

    float mul;

    lbl1:
    {
        while (j < (colms-1))
        { 
            for (k=0,l=j; k<(rows) && l<(colms-1); k++,l++)
            {
                if ( rref[k+n][l] != 1 )
                {
                    j++;
                    if ((k+n) == (rows-1)) n++;
                    goto lbl1;
                }

                else
                {
                    mul = rref[k][l];
                    for (m=l; m<colms; m++)
                        rref[k][m] = rref[k][m] -(mul*rref[k+n][m]);
                }
            }
        j++;
        n++;
        }
    }
}

unsigned int lSearch(char *str, char ch)
{
    unsigned int freq = 0;

    while (*str != '\0')
    {
        if (*str == ch)
            freq++;

        str = str + 1;
    }

    return freq;
}