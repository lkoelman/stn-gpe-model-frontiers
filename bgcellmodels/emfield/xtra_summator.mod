COMMENT

Sum LFP sources from other compartments and assign the summed LFP signal
to a recordable variable.

Written by Lucas Koelman

ENDCOMMENT


NEURON {
    POINT_PROCESS xtra_sum
    POINTER temp_ptr
    POINTER donotuse_sources
    RANGE summed
}


VERBATIM

// Linked list node for storing refs to observed hoc variables
typedef struct node {
    double* hoc_ref;    // hoc reference to observed LFP source variable
    struct node* next;  // next node in linked list
} LfpSource;

ENDVERBATIM


UNITS {
    (nV) = (nanovolt)
}



ASSIGNED {
    temp_ptr
    donotuse_sources
    summed (nV)
}


CONSTRUCTOR {
VERBATIM {
    // Snippet based on pattern.mod
    LfpSource** sources = (LfpSource**)(&(_p_donotuse_sources));

    LfpSource* first_node = emalloc(sizeof(LfpSource));
    first_node->hoc_ref = NULL;
    first_node->next = NULL;
    
    *sources = first_node;

}
ENDVERBATIM
}


DESTRUCTOR {
VERBATIM {

    // Free memory occupied by linked list
    LfpSource* this_node = (LfpSource*)(_p_donotuse_sources);
    LfpSource* next_node;
    while (this_node != NULL) {
        next_node = this_node->next;
        free(this_node);
        this_node = next_node;
    }

}
ENDVERBATIM
}


: Add observed LFP source.
:
: PYTHON USAGE
: ------------
: 
: >>> summator = LfpSummator(soma(0.5))
: >>> for sec in h.allsec():
: >>>     for seg in sec:
: >>>         h.setpointer(seg._ref_i_membrane, 'temp_ptr', summator)
: >>>         summator.add_source()
:
:
FUNCTION add_imemb_source() {
VERBATIM
    
    // Look for end of linked list and append observed variable
    LfpSource* current = (LfpSource*)(_p_donotuse_sources);
    while (current->hoc_ref != NULL) {
        current = current->next;
    }
    current->hoc_ref = _p_temp_ptr;

    // Allocate node for next call
    current->next = emalloc(sizeof(LfpSource)); // for next call
    current->next->hoc_ref = NULL;
    current->next->next = NULL;

    // fprintf(stderr, "Added ref to group %d\n", group_id);
ENDVERBATIM
}


: Update pointer of imemb source, e.g. after it has become invalid.
FUNCTION update_imemb_ptr(index) {
VERBATIM
    
    LfpSource* current = (LfpSource*)(_p_donotuse_sources);
    unsigned int current_index = 0;
    while (current_index < _lindex) {
        current = current->next;
        ++current_index;
    }
    if (current == NULL) {
        hoc_execerror("Index too high: no pointer at index", 0);
    }
    current->hoc_ref = _p_temp_ptr;

ENDVERBATIM
}


: Sum all observed LFP sources and assign to 'summed'
: NOTE: AFTER SOLVE is called multiple times per step -> leads to 42% increase
:       in simulation time in test
BREAKPOINT {
VERBATIM
    LfpSource* current = (LfpSource*)(_p_donotuse_sources);
    summed = 0.0;
    while (current != NULL) {
        if (current->hoc_ref != NULL) {
            double i_membrane = *(current->hoc_ref);
            summed += i_membrane;
        }
        current = current->next;
    }
ENDVERBATIM
}