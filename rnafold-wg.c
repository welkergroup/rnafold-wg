/*
 * Compile, az rnafold examples mappában (make után):
 *     gcc -Wall -fopenmp -I ../src -o rnafold-wg rnafold-wg.c -L ../src/ViennaRNA/.libs -lRNA -lm
 */

#define _GNU_SOURCE

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/params/basic.h>
#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/eval.h>
#include <ViennaRNA/fold.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/subopt.h>

#define SPACER_LENGTH 20
#define RESULT_ROW_LENGTH (100 + SPACER_LENGTH)

FLT_OR_DBL kT = 0.616321;
int radius = 8;
char *rar_constraint = "....................((((((..((((....))))....))))))....................................................";
char *sl1_constraint = "....................................................(((..).)).........................................";
char *sl2_constraint = "....................................................................((((....))))......................";
char *sl3_constraint = ".................................................................................((((((...))))))......";
char *mml_postfix = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTT";
char *ngg_postfix = "GGGTTAGAGCTAGAAATAGCAAGTTAACCTAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTT";

typedef struct linked_list_s linked_list_t;

struct linked_list_s {
    linked_list_t * next;
    char * spacer;
    char * result_row;
};

typedef struct subopt_state_s {
    FLT_OR_DBL free_energy_of_ensemble;
    FLT_OR_DBL summed_percent;
} subopt_state_t;

void
subopt_callback(const char *structure,
                float energy,
                void *data) {
    if (structure) {
        subopt_state_t *state = (subopt_state_t *) data;
        FLT_OR_DBL percent = 100.0 * exp((state->free_energy_of_ensemble - energy) / kT);
        state->summed_percent += percent;
    }
}

FLT_OR_DBL
calculate_sequence(const char *sequence, const char *constraint) {
    char *mfe_structure = vrna_alloc(sizeof(char) * (strlen(sequence) + 1));
    char *prob_string = vrna_alloc(sizeof(char) * (strlen(sequence) + 1));

    vrna_md_t md;

    vrna_md_set_default(&md);
    md.dangles = 2;
    md.noLP = 1;
    md.compute_bpp = 0;
    md.uniq_ML = 1;
    md.min_loop_size = 2;

    /* get a vrna_fold_compound with default settings */
    vrna_fold_compound_t *vc = vrna_fold_compound(sequence, &md, VRNA_OPTION_DEFAULT);

    /* call MFE function */
    double mfe = (double) vrna_mfe(vc, mfe_structure);

    /* rescale parameters for Boltzmann factors */
    vrna_exp_params_rescale(vc, &mfe);

    /* call PF function */
    FLT_OR_DBL en = vrna_pf(vc, prob_string);

    subopt_state_t state = {
            // .free_energy_of_ensemble = roundf(en * 100) / 100,
            .free_energy_of_ensemble = en,
            .summed_percent =  0
    };

    if (constraint != NULL) {
        vrna_constraints_add(vc, constraint, VRNA_CONSTRAINT_DB_DEFAULT);
    }

    vrna_subopt_cb(vc, radius * 100, &subopt_callback, (void *) &state);

    /* free pseudo dot-bracket probability string */
    free(prob_string);

    /* free mfe structure */
    free(mfe_structure);

    /* free memory occupied by vrna_fold_compound */
    vrna_fold_compound_free(vc);

    return state.summed_percent;
}

char * calculate_spacer(const char *spacer, const char *postfix, const char *constraint) {
    int length = strlen(spacer);
    if (length != SPACER_LENGTH) {
        printf("invalid length of %d for %s\n", length, spacer);
        return NULL;
    }

    char *sequence = vrna_alloc(sizeof(char) * (length + strlen(postfix) + 1));
    strcpy(sequence, spacer);
    strcpy(sequence + length, postfix);

    if (constraint == NULL) {
        FLT_OR_DBL full_percent = calculate_sequence(sequence, NULL);
        FLT_OR_DBL rar_percent = calculate_sequence(sequence, rar_constraint);
        FLT_OR_DBL sl1_percent = calculate_sequence(sequence, sl1_constraint);
        FLT_OR_DBL sl2_percent = calculate_sequence(sequence, sl2_constraint);
        FLT_OR_DBL sl3_percent = calculate_sequence(sequence, sl3_constraint);

        free(sequence);

        char *result_row = vrna_alloc(sizeof(char) * RESULT_ROW_LENGTH);
        snprintf(result_row, RESULT_ROW_LENGTH, "%s,%.12f,%.12f,%.12f,%.12f,%.12f\n", spacer, full_percent, rar_percent, sl1_percent, sl2_percent, sl3_percent);
        return result_row;
    } else {
        FLT_OR_DBL full_percent = calculate_sequence(sequence, NULL);
        FLT_OR_DBL constraint_percent = calculate_sequence(sequence, constraint);

        free(sequence);

        char *result_row = vrna_alloc(sizeof(char) * RESULT_ROW_LENGTH);
        snprintf(result_row, RESULT_ROW_LENGTH, "%s,%.12f,%.12f\n", spacer, full_percent, constraint_percent);
        return result_row;
    }
}

int
main(int argc, char *argv[]) {
    const char *postfix = NULL;
    const char *constraint = NULL;

    if (argc == 1) {
        postfix = mml_postfix;
    }
    else if (argc == 2) {
        if (strcmp(argv[1], "wt") == 0) {
            postfix = mml_postfix;
        }
        else if (strcmp(argv[1], "ngg") == 0) {
            postfix = ngg_postfix;
        }
        else if (strlen(argv[1]) == 82) {
            postfix = argv[1];
        }
        else {
            fprintf(stderr, "No guide sequence specified -- no arguments or guide sequence expected (wt, ngg, or 82 chars)\n");
            return EXIT_FAILURE;
        }
    }
    else if (argc == 3) {
        postfix = argv[1];
        constraint = argv[2];

        if (strlen(postfix) + SPACER_LENGTH != strlen(constraint)) {
            fprintf(stderr, "Postfix and constraint length mistmatch (%li + %i != %li)\n", strlen(postfix), SPACER_LENGTH, strlen(constraint));
            return EXIT_FAILURE;
        }
    } else {
            fprintf(stderr, "Usage: %s [guide [constraint]]\n", argv[0]);
            return EXIT_FAILURE;
    }

    if (constraint == NULL) {
        fprintf(stderr, "Using %s guide, reading from stdin.\n", postfix == ngg_postfix ? "ngg" : postfix == mml_postfix ? "wt" : "custom");
    }
    else {
        fprintf(stderr, "Using custom guide with custom constraint, reading from stdin.\n");
    }

    linked_list_t *results = NULL;

#pragma omp parallel
#pragma omp single
    {
        linked_list_t *tail = results;

        char buf[23];
        while (fgets(buf, 22, stdin) != 0) {
            buf[strlen(buf) - 1] = '\0';

            char *spacer = vrna_alloc(sizeof(char) * (strlen(buf) + 1));
            strcpy(spacer, buf);

            linked_list_t *ll = vrna_alloc(sizeof(linked_list_t)); 
            ll->next = NULL;
            ll->spacer = spacer;

            if (tail != NULL) {
                tail->next = ll;
                tail = ll;
            } else {
                tail = results = ll;
            }

#pragma omp task firstprivate(ll, spacer) shared(postfix)
            {
                int threadId = 0;
                char name[16] = { 0 };

                #ifdef _OPENMP
                threadId = omp_get_thread_num();
                #endif

                snprintf(name, sizeof(name) - 1, "omp-%03d", threadId);

                pthread_setname_np(pthread_self(), name);

                ll->result_row = calculate_spacer(spacer, postfix, constraint);
            }
        }
    }
#pragma omp taskwait

    if (constraint == NULL) {
        printf("spacer,full_percent,rar_percent,sl1_percent,sl2_percent,sl3_percent\n");
    } else {
        printf("spacer,full_percent,constraint_percent\n");
    }

    linked_list_t *next = results;
    while (next != NULL) {
        linked_list_t *cur = next;
        next = cur->next;

        if (cur->result_row != NULL) {
            printf("%s", cur->result_row);
            free(cur->result_row);
        } else {
            printf("%s,,,,,\n", cur->spacer);
        }
        free(cur->spacer);
        free(cur);
    }

    return EXIT_SUCCESS;
}
