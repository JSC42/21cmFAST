
void CopyArray(double *x1, double *x0)
{
    // Copy contents in x0 to x1, both must be of the same size
    int idx;
    for (idx=0; idx < Redshift_Size*Ek_axis_Size; idx++)
    {
        x1[idx] = x0[idx];
    }
}
