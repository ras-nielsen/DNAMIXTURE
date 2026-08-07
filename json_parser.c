#include "json_parser.h"

// Global variable
PopFreqHashTable *global_popfreq_table = NULL;

int string_to_nucleotide(const char *str)
{
        // Skip leading whitespace
        while (*str && (*str == ' ' || *str == '\t' || *str == '\n' || *str == '\r')) {
                str++;
        }

        if (!*str) {
                fprintf(stderr, "Error: Empty nucleotide string\n");
                exit(1);
        }

        return nucleotide_to_int(*str);
}

// Helper function: load entire file into memory
char* load_file(const char *filename, size_t *file_size)
{
        FILE *f = fopen(filename, "rb");
        if (!f) {
                fprintf(stderr, "Error: Cannot open file '%s'\n", filename);
                perror("fopen");
                exit(1);
        }

        // Get file size
        fseek(f, 0, SEEK_END);
        *file_size = ftell(f);
        fseek(f, 0, SEEK_SET);

        // Allocate buffer
        char *buffer = (char*)malloc(*file_size + 1);
        if (!buffer) {
                fprintf(stderr, "Error: Cannot allocate memory for file (size: %zu bytes)\n", *file_size);
                exit(1);
        }

        // Read file
        size_t read_size = fread(buffer, 1, *file_size, f);
        if (read_size != *file_size) {
                fprintf(stderr, "Error: Failed to read entire file\n");
                exit(1);
        }
        buffer[*file_size] = '\0';  // Null terminate

        fclose(f);
        return buffer;
}

// Helper function: find substring in buffer
char* find_string(char *haystack, const char *needle)
{
        return strstr(haystack, needle);
}

// Helper function: skip whitespace
char* skip_whitespace(char *p)
{
        while (*p && (*p == ' ' || *p == '\t' || *p == '\n' || *p == '\r')) {
                p++;
        }
        return p;
}

// Helper function: parse a quoted string, returns newly allocated string
char* parse_quoted_string(char **p)
{
        *p = skip_whitespace(*p);
        if (**p != '"') {
                fprintf(stderr, "Error: Expected '\"' but found '%c' (ASCII %d)\n", **p, (int)**p);
                fprintf(stderr, "Context: %.50s\n", *p);
                exit(1);
        }
        (*p)++;  // Skip opening quote

        char *start = *p;
        while (**p && **p != '"') {
                if (**p == '\\') {
                        (*p)++;  // Skip escape char
                        if (**p) (*p)++;  // Skip escaped char
                } else {
                        (*p)++;
                }
        }

        if (**p != '"') {
                fprintf(stderr, "Error: Unterminated string\n");
                exit(1);
        }

        size_t len = *p - start;
        char *result = (char*)malloc(len + 1);
        strncpy(result, start, len);
        result[len] = '\0';

        (*p)++;  // Skip closing quote
        return result;
}

// Helper function: parse a number (integer or float)
double parse_number(char **p)
{
        *p = skip_whitespace(*p);
        char *end;
        double value = strtod(*p, &end);
        if (end == *p) {
                fprintf(stderr, "Error: Expected number\n");
                exit(1);
        }
        *p = end;
        return value;
}

// Hash function for position keys
unsigned int hash_position_key(const char *key, int table_size)
{
        unsigned int hash = 5381;
        int c;
        while ((c = *key++)) {
                hash = ((hash << 5) + hash) + c; // hash * 33 + c
        }
        return hash % table_size;
}

// Create hash table
PopFreqHashTable* create_popfreq_table(int size)
{
        PopFreqHashTable *table = (PopFreqHashTable*)malloc(sizeof(PopFreqHashTable));
        table->size = size;
        table->count = 0;
        table->buckets = (PopFreqNode**)calloc(size, sizeof(PopFreqNode*));
        return table;
}

// Insert into hash table
void insert_popfreq(PopFreqHashTable *table, const char *position_key, double AF[4])
{
        unsigned int index = hash_position_key(position_key, table->size);

        PopFreqNode *node = (PopFreqNode*)malloc(sizeof(PopFreqNode));
        node->position_key = strdup(position_key);
        memcpy(node->AF, AF, 4 * sizeof(double));
        node->next = table->buckets[index];
        table->buckets[index] = node;
        table->count++;
}

// Lookup in hash table
int lookup_popfreq(PopFreqHashTable *table, const char *position_key, double AF[4])
{
        if (!table) return 0;

        unsigned int index = hash_position_key(position_key, table->size);
        PopFreqNode *node = table->buckets[index];

        while (node) {
                if (strcmp(node->position_key, position_key) == 0) {
                        memcpy(AF, node->AF, 4 * sizeof(double));
                        return 1;
                }
                node = node->next;
        }
        return 0;
}

// Free hash table
void free_popfreq_table(PopFreqHashTable *table)
{
        if (!table) return;

        int i;
        for (i = 0; i < table->size; i++) {
                PopFreqNode *node = table->buckets[i];
                while (node) {
                        PopFreqNode *next = node->next;
                        free(node->position_key);
                        free(node);
                        node = next;
                }
        }
        free(table->buckets);
        free(table);
}

// Parse population_freqs and store in hash table
// Structure: {"10:95662": {"G": 0.65, "C": 0.35}, ...}
void parse_population_freqs(char *json_start)
{
        fprintf(stderr, "Parsing population_freqs into hash table...\n");

        char *p = find_string(json_start, "\"population_freqs\"");
        if (!p) {
                fprintf(stderr, "Warning: Cannot find 'population_freqs' section\n");
                global_popfreq_table = NULL;
                return;
        }

        // Find the opening brace
        p = strchr(p, '{');
        if (!p) {
                fprintf(stderr, "Error: Malformed population_freqs section\n");
                exit(1);
        }
        p++;  // Skip '{'

        // Create hash table (size ~4M for 2.8M entries, load factor ~0.7)
        global_popfreq_table = create_popfreq_table(4000000);

        int parsed_count = 0;

        // Parse each position entry
        while (1) {
                p = skip_whitespace(p);
                if (*p == '}') break;  // End of population_freqs
                if (*p == ',') {
                        p++;
                        continue;
                }

                // Parse position key like "10:95662"
                if (*p != '"') break;
                char *position_key = parse_quoted_string(&p);

                // Skip colon separator
                p = skip_whitespace(p);
                if (*p != ':') {
                        free(position_key);
                        break;
                }
                p++;

                // Parse frequency object {"G": 0.65, "C": 0.35}
                p = skip_whitespace(p);
                if (*p != '{') {
                        free(position_key);
                        break;
                }
                p++;

                // Parse allele frequencies
                double AF[4];
                int i;
                for (i = 0; i < 4; i++) AF[i] = 0.001;  // Initialize

                while (1) {
                        p = skip_whitespace(p);
                        if (*p == '}') {
                                p++;  // Skip closing brace
                                break;
                        }
                        if (*p == ',') {
                                p++;
                                continue;
                        }

                        // Parse nucleotide key
                        if (*p != '"') break;
                        char *nuc_str = parse_quoted_string(&p);

                        // Skip colon
                        p = skip_whitespace(p);
                        if (*p == ':') p++;

                        // Parse frequency value
                        double freq = parse_number(&p);

                        // Store in AF array
                        int nuc_idx = nucleotide_to_int(nuc_str[0]);
                        if (nuc_idx >= 0 && nuc_idx < 4) {
                                AF[nuc_idx] = freq;
                        }
                        free(nuc_str);
                }

                // Renormalize
                double sum = 0.0;
                for (i = 0; i < 4; i++) sum += AF[i];
                for (i = 0; i < 4; i++) AF[i] /= sum;

                // Insert into hash table
                insert_popfreq(global_popfreq_table, position_key, AF);
                free(position_key);

                parsed_count++;
                if (parsed_count % 100000 == 0) {
                        fprintf(stderr, "  Parsed %d population frequencies...\r", parsed_count);
                }
        }

        fprintf(stderr, "\nSuccessfully parsed %d population frequencies into hash table\n", parsed_count);
}

// Look up allele frequencies for a specific position from hash table
// Returns 1 on success, 0 if position not found
int lookup_population_freqs(char *json_start, const char *position_key, double AF[4])
{
        // Simply look up in hash table - O(1) average case
        return lookup_popfreq(global_popfreq_table, position_key, AF);
}

// Count number of positions in position_reads by counting top-level entries
int count_positions(char *json_start)
{
        char *p = find_string(json_start, "\"position_reads\"");
        if (!p) {
                return 0;
        }

        p = strchr(p, '{');
        if (!p) return 0;
        p++;  // Skip opening '{'

        // Count comma-separated entries at depth 1
        int count = 0;
        int brace_depth = 0;
        int in_string = 0;
        int found_first = 0;

        while (*p) {
                if (*p == '"' && (p == json_start || *(p-1) != '\\')) {
                        in_string = !in_string;
                } else if (!in_string) {
                        if (*p == '{') {
                                brace_depth++;
                        } else if (*p == '}') {
                                if (brace_depth == 0) break;  // End of position_reads
                                brace_depth--;
                        } else if (*p == ',' && brace_depth == 0) {
                                count++;
                        } else if (brace_depth == 0 && !found_first && *p != ' ' && *p != '\n' && *p != '\r' && *p != '\t') {
                                found_first = 1;
                        }
                }
                p++;
        }

        if (found_first) count++;  // Count the last entry

        return count;
}

// Parse a single position entry from position_reads
// Returns pointer after this entry, or NULL on error
char* parse_position_entry(char *p, int snp_idx, char *json_start)
{
        // We're at the start of a key like "10:1978952"
        p = skip_whitespace(p);
        if (*p != '"') {
                return NULL;
        }

        char *key_start = parse_quoted_string(&p);

        // Save the original position key for population_freqs lookup
        char *position_key = strdup(key_start);

        // Parse chrom:pos
        char *colon = strchr(key_start, ':');
        if (!colon) {
                free(key_start);
                free(position_key);
                return NULL;
        }

        *colon = '\0';
        char *chrom = key_start;
        int pos = atoi(colon + 1);

        // Skip the ':' between key and value
        p = skip_whitespace(p);
        if (*p != ':') {
                free(key_start);
                free(position_key);
                return NULL;
        }
        p++;

        // Now parse the value object {"reads": [...], "suspect_gt": [...], "victim_gt": [...]}
        p = skip_whitespace(p);
        if (*p != '{') {
                free(key_start);
                free(position_key);
                return NULL;
        }
        p++;  // Skip '{'

        // Initialize arrays for this SNP
        int *reads_nt = NULL;
        double *reads_qual = NULL;
        int num_reads = 0;
        int suspect_gt[2] = {-1, -1};
        int victim_gt[2] = {-1, -1};

        // Parse the object fields
        while (1) {
                p = skip_whitespace(p);
                if (*p == '}') {
                        p++;  // Skip closing '}'
                        break;
                }
                if (*p == ',') {
                        p++;
                        continue;
                }

                // Parse field name
                char *field_name = parse_quoted_string(&p);

                p = skip_whitespace(p);
                if (*p != ':') {
                        free(field_name);
                        free(key_start);
                        free(position_key);
                        return NULL;
                }
                p++;  // Skip ':'

                p = skip_whitespace(p);

                if (strcmp(field_name, "reads") == 0) {
                        // Parse reads array: [[" G", 37], ["A", 37], ...]
                        if (*p != '[') {
                                free(field_name);
                                free(key_start);
                                free(position_key);
                                return NULL;
                        }
                        p++;  // Skip opening '[' of reads array

                        // Count reads first
                        // We've already skipped the opening '[', so we're inside the array
                        // Each read is a nested array like ["G", 37]
                        // Count opening '[' brackets until we hit the closing ']' of the reads array
                        char *p_count = p;
                        num_reads = 0;
                        int bracket_depth = 0;  // We're inside the outer array, so depth 0 = inside reads array
                        while (*p_count) {
                                if (*p_count == '[') {
                                        bracket_depth++;
                                        if (bracket_depth == 1) {
                                                // This is an opening '[' of a read pair like ["G", 37]
                                                num_reads++;
                                        }
                                } else if (*p_count == ']') {
                                        bracket_depth--;
                                        if (bracket_depth < 0) {
                                                // This ']' closes the reads array itself
                                                break;
                                        }
                                }
                                p_count++;
                        }

                        // Allocate arrays
                        reads_nt = (int*)malloc(num_reads * sizeof(int));
                        reads_qual = (double*)malloc(num_reads * sizeof(double));

                        // Parse each read
                        int read_idx = 0;
                        while (read_idx < num_reads) {
                                p = skip_whitespace(p);
                                if (*p == ']') break;  // End of reads array
                                if (*p == ',') {
                                        p++;
                                        continue;
                                }

                                if (*p != '[') break;
                                p++;  // Skip opening '[' of read pair

                                // Parse nucleotide
                                char *nt_str = parse_quoted_string(&p);
                                reads_nt[read_idx] = string_to_nucleotide(nt_str);
                                free(nt_str);

                                p = skip_whitespace(p);
                                if (*p == ',') p++;

                                // Parse quality
                                reads_qual[read_idx] = parse_number(&p);

                                p = skip_whitespace(p);
                                if (*p == ']') p++;  // Skip closing ']' of read pair

                                read_idx++;
                        }

                        // Now we should be at the closing ']' of the entire reads array
                        p = skip_whitespace(p);
                        if (*p == ']') p++;  // Skip closing ']' of reads array

                } else if (strcmp(field_name, "suspect_gt") == 0 || strcmp(field_name, "victim_gt") == 0) {
                        // Parse genotype array: ["G", "G"]
                        if (*p != '[') {
                                free(field_name);
                                free(key_start);
                                free(position_key);
                                if (reads_nt) free(reads_nt);
                                if (reads_qual) free(reads_qual);
                                return NULL;
                        }
                        p++;  // Skip '['

                        char *gt1_str = parse_quoted_string(&p);
                        p = skip_whitespace(p);
                        if (*p == ',') p++;
                        char *gt2_str = parse_quoted_string(&p);

                        int gt1 = string_to_nucleotide(gt1_str);
                        int gt2 = string_to_nucleotide(gt2_str);
                        free(gt1_str);
                        free(gt2_str);

                        if (strcmp(field_name, "suspect_gt") == 0) {
                                suspect_gt[0] = gt1;
                                suspect_gt[1] = gt2;
                        } else {
                                victim_gt[0] = gt1;
                                victim_gt[1] = gt2;
                        }

                        p = skip_whitespace(p);
                        if (*p == ']') p++;  // Skip closing ']'
                }

                free(field_name);
        }

        // Now populate the SNPData structure
        if (num_reads > 0 && suspect_gt[0] >= 0 && victim_gt[0] >= 0) {
                allocate_snp_reads(&global_snp_data->snps[snp_idx], num_reads);
                set_snp_genotypes(&global_snp_data->snps[snp_idx],
                                 victim_gt[0], victim_gt[1], suspect_gt[0], suspect_gt[1]);

                // Set reads
                int i;
                for (i = 0; i < num_reads; i++) {
                        set_snp_read(&global_snp_data->snps[snp_idx], i, reads_nt[i], reads_qual[i]);
                }

                // Set allele frequencies - look up from population_freqs
                double af[4];
                if (!lookup_population_freqs(json_start, position_key, af)) {
                        // If position not found in population_freqs, use uniform frequencies
                        fprintf(stderr, "Warning: Position %s not found in population_freqs, using uniform frequencies\n", position_key);
                        af[0] = af[1] = af[2] = af[3] = 0.25;
                }
                set_snp_af(&global_snp_data->snps[snp_idx], af);
        }

        // Cleanup
        if (reads_nt) free(reads_nt);
        if (reads_qual) free(reads_qual);
        free(key_start);
        free(position_key);

        return p;
}

int parse_data(){
        fprintf(stderr, "Loading JSON file: %s\n", global_opts.infile);

        size_t file_size;
        char *json = load_file(global_opts.infile, &file_size);

        fprintf(stderr, "File loaded: %zu bytes\n", file_size);

        // Parse population_freqs
        parse_population_freqs(json);

        // Count positions in position_reads
        int num_positions = count_positions(json);
        fprintf(stderr, "Found %d positions in position_reads\n", num_positions);

        if (num_positions == 0) {
                fprintf(stderr, "Error: No positions found in position_reads\n");
                exit(1);
        }

        // Allocate MultiSNPData
        global_snp_data = allocate_multisnp_data(num_positions);

        fprintf(stderr, "Parsing position_reads section...\n");

        // Find position_reads section
        char *p = find_string(json, "\"position_reads\"");
        if (!p) {
                fprintf(stderr, "Error: Cannot find position_reads section\n");
                exit(1);
        }

        p = strchr(p, '{');
        if (!p) {
                fprintf(stderr, "Error: Malformed position_reads section\n");
                exit(1);
        }
        p++;  // Skip '{'

        // Parse each position entry
        int snp_idx = 0;
        int parsed_count = 0;
        while (snp_idx < num_positions && *p) {
                p = skip_whitespace(p);
                if (*p == '}') break;  // End of position_reads
                if (*p == ',') {
                        p++;
                        continue;
                }

                char *next_p = parse_position_entry(p, snp_idx, json);
                if (next_p == NULL) {
                        fprintf(stderr, "Error parsing position %d\n", snp_idx);
                        break;
                }
                p = next_p;
                parsed_count++;
                snp_idx++;

                if (parsed_count % 500 == 0) {
                        fprintf(stderr, "  Parsed %d / %d positions...\r", parsed_count, num_positions);
                }
        }

        fprintf(stderr, "\nSuccessfully parsed %d positions\n", parsed_count);

        free(json);
        return 0;
}
