#ifndef COMMS_H
#define COMMS_H
#include <mpi.h>

/*
Class to define MPI communications between with a neighbouring processor
*/
class Comms
{
public:
    int neighbour_proc;
    int comm_inds[2][4];        // contains (start_i, start_j, end_i, end_j) for sends [0][:] and recieves [1][:]
    MPI_Datatype comm_types[2]; // send_type [0] and recv_type [1]

    // Create custom MPI_Datatypes for sending and receiving boundary grid data from a neighbouring processor
    void create_MPI_types(bool **proc_grid);
    // Explicitly free existing MPI Datatype's
    void free_types();

private:
    bool comm_with_self = true; // gets set to false if create_MPI_type() is called
};

#endif