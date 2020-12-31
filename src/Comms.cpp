#include "../include/Comms.h"
#include <vector>

/*
Create custom MPI_Datatypes for sending and receiving boundary grid data from a neighbouring processor
*/
void Comms::create_MPI_types(bool **proc_grid)
{
    std::vector<int> block_lengths;
    std::vector<MPI_Aint> addresses;
    std::vector<MPI_Datatype> typelist;

    for (int n = 0; n < 2; n++)
    {
        MPI_Datatype MPI_type;

        int start_i = comm_inds[n][0];
        int start_j = comm_inds[n][1];
        int end_i = comm_inds[n][2];
        int end_j = comm_inds[n][3];

        int j_length = end_j - start_j + 1;

        for (int i = start_i; i <= end_i; i++)
        {
            block_lengths.push_back(j_length);
            MPI_Aint temp;
            MPI_Get_address(&proc_grid[i][start_j], &temp);
            addresses.push_back(temp);
            typelist.push_back(MPI_C_BOOL);
        }

        MPI_Type_create_struct(typelist.size(), block_lengths.data(), addresses.data(), typelist.data(), &MPI_type);
        MPI_Type_commit(&MPI_type);
        comm_types[n] = MPI_type;

        block_lengths.clear();
        addresses.clear();
        typelist.clear();
    }
    comm_with_self = false;
}

/*
Explicitly free existing MPI Datatype's
*/
void Comms::free_types()
{
    if (!comm_with_self)
    {
        MPI_Type_free(&comm_types[0]);
        MPI_Type_free(&comm_types[1]);
    }
}