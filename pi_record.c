#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>
#include <locale.h>
#include <string.h>

char c;
#define clear while((c = getchar()) != '\n') {}
#define mpfr_null (mpfr_ptr)0
#define prec_buff 1000
// métodos que calculam uma aproximação de pi usando as
// fórmulas de K. Takano e F. C. M. Stormer, respectivamente,
// com precisão {prec} e guardam o resultado em {out} e em um arquivo
int K_Takano_method(mpfr_t out);
int FCM_Stormer_method(mpfr_t out);
// gera uma string de formatação bonitinha para o mpfr_printf que
// inclui uma precisão variável
char *generate_format_str(char *label, long int prec);
// compara nossa estimativa A com a estimativa B e retorna o número de
// casas decimais equivalentes entre as duas estimativas (nossa precisão)
long int compare_results(mpfr_t Takano_result, mpfr_t Stormer_result);

// nossa precisão base (a verdadeira precisão é menor por conta das limitações das fórmulas)
long int prec = -1;

int main(int argc, char const *argv[]) {
	setlocale(LC_ALL, "");

	if(argc > 1) {
		prec = atol(argv[1]);
	} else {
		while(prec == -1) {
			printf("com qual precisão você quer calcular pi (em número de casas decimais)?\n");
			scanf("%ld", &prec);
			clear;
		}
	}
	// só pra garantir maior corretude quando a precisão fornecida pelo usuário é baixa de mais
	prec += prec_buff;

	// usando nosso método para aproximar pi com a fórmula de K. Takano
	mpfr_t pi_K_Takano;
	mpfr_init2(pi_K_Takano, prec);
	K_Takano_method(pi_K_Takano);
	printf("\n");
	// e agora com a fórmula de F. C. M. Stormer
	mpfr_t pi_FCM_Stormer;
	mpfr_init2(pi_FCM_Stormer, prec);
	FCM_Stormer_method(pi_FCM_Stormer);

	// colocando o valor de pi do K. Takano em um arquivo, só por curiosidade
	// e para termos um registro "permanente", diferentemente do que é printado no terminal
	FILE *output_file = fopen("pi_K_Takano.txt", "w");
	if (output_file == NULL) {
	    perror("Erro ao criar arquivo de output");
	    return 1;
	}
	mpfr_out_str(output_file, 10, prec - prec_buff + 1, pi_K_Takano, MPFR_RNDZ); // o - prec_buff é pra compensar a precisão extra da linha 40
	fclose(output_file);

	// colocando o valor de pi do F. C. M. Stormer em outro arquivo
	output_file = fopen("pi_FCM_Stormer.txt", "w");
	if (output_file == NULL) {
	    perror("Erro ao criar arquivo de output");
	    return 1;
	}
	mpfr_out_str(output_file, 10, prec - prec_buff + 1, pi_FCM_Stormer, MPFR_RNDZ);
	fclose(output_file);

	// achando quantas casas decimais estão equivalentes nas duas estimativas e printando elas
	long int correct_decimal_places = compare_results(pi_K_Takano, pi_FCM_Stormer);
	char pi_correct_portion[correct_decimal_places + 5];
	pi_correct_portion[0] = '\0';
	mpfr_sprintf(pi_correct_portion, generate_format_str("", correct_decimal_places + 2), pi_K_Takano); // o +2 impede o generate_format_str de arredondar a parte que queremos
	pi_correct_portion[strlen(pi_correct_portion) - 2] = '\0';
	printf("\nPorção equivalente das aproximações: %s\n\n", pi_correct_portion);
	printf("RESULTADO:\n%ld casas decimais estimadas\n%ld casas decimais equivalentes\n", prec - prec_buff, correct_decimal_places);

	return 0;
}

int K_Takano_method(mpfr_t out) {
	// criando listas de floats, racionais e ints de precisão arbitrária
	// para os valores que serão usados na fórmula
	mpfr_t arctg_f[4];
	mpfr_t frac_f[4];
	mpfr_t term_f[4];
	mpq_t frac_q[4];
	mpz_t coef_z[4];
	int denominators[4] = {49, 57, 239, 110443};
	int coefficients[4] = {12 * 4, 32 * 4, (-5) * 4, 12 * 4};

	for(int i = 0; i < 4; i++) {
		// inicializando as variáveis
		mpfr_init2(arctg_f[i], prec);
		mpfr_init2(frac_f[i], prec);
		mpfr_init2(term_f[i], prec);
		mpq_init(frac_q[i]);
		mpz_init(coef_z[i]);
		// setando as frações
		mpq_set_ui(frac_q[i], 1, denominators[i]);
		// convertendo elas para mpfr_t
		mpfr_set_q(frac_f[i], frac_q[i], MPFR_RNDN);
		// setando os valores das arco-tangentes
		mpfr_atan(arctg_f[i], frac_f[i], MPFR_RNDN);
		// setando os valores dos coeficientes
		mpz_set_si(coef_z[i], coefficients[i]);
		// setando os termos inicialmente como apenas os coeficientes
		mpfr_set_z(term_f[i], coef_z[i], MPFR_RNDN);
		// multiplicando nossos coeficientes pelas arco-tangentes
		mpfr_mul(term_f[i], term_f[i], arctg_f[i], MPFR_RNDN);

	}

	// colocando o valor do primeiro termo no output
	mpfr_swap(out, term_f[0]);
	// somando/subtraindo os outros termos
	for(int i = 1; i < 4; i++) mpfr_add(out, out, term_f[i], MPFR_RNDN);

	//gerando uma string de formatação personalizada para o mpfr_printf
	char *format_str = generate_format_str("PI (pelo método de K. Takano): ", prec - prec_buff + 2);
	// printando nosso resultado
	char Takano_result[prec + 5];
	mpfr_sprintf(Takano_result, format_str, out);
	//Takano_result[strlen(Takano_result) - 3] = '\n';
	Takano_result[strlen(Takano_result) - 2] = '\0';
	printf("%s", Takano_result);
	printf("\n");
	// liberando a string de formatação da memória
	free(format_str);

	for(int i = 0; i < 4; i++) {
		// limpando as variáveis usadas
		mpfr_clear(arctg_f[i]);
		mpfr_clear(term_f[i]);
		mpfr_clear(frac_f[i]);
		mpz_clear(coef_z[i]);
		mpq_clear(frac_q[i]);
	}
	return 0;
}

int FCM_Stormer_method(mpfr_t out) {
	// tudo aqui funciona da mesma forma que na outra função
	// só o que muda são os coeficientes e as frações
	mpfr_t arctg_f[4];
	mpfr_t frac_f[4];
	mpfr_t term_f[4];
	mpq_t frac_q[4];
	mpz_t coef_z[4];
	int denominators[4] = {57, 239, 682, 12943};
	int coefficients[4] = {44 * 4, 7 * 4, (-12) * 4, 24 * 4};

	for(int i = 0; i < 4; i++) {
		mpfr_init2(arctg_f[i], prec);
		mpfr_init2(frac_f[i], prec);
		mpfr_init2(term_f[i], prec);
		mpq_init(frac_q[i]);
		mpz_init(coef_z[i]);

		mpq_set_ui(frac_q[i], 1, denominators[i]);
		mpfr_set_q(frac_f[i], frac_q[i], MPFR_RNDN);
		mpfr_atan(arctg_f[i], frac_f[i], MPFR_RNDN);
		mpz_set_si(coef_z[i], coefficients[i]);
		mpfr_set_z(term_f[i], coef_z[i], MPFR_RNDN);
		mpfr_mul(term_f[i], term_f[i], arctg_f[i], MPFR_RNDN);
	}

	mpfr_swap(out, term_f[0]);
	for(int i = 1; i < 4; i++) mpfr_add(out, out, term_f[i], MPFR_RNDN);

	char Stormer_result[prec + 5];
	char *format_str = generate_format_str("PI (pelo método de F. C. M. Stormer): ", prec - prec_buff + 2);
	mpfr_sprintf(Stormer_result, format_str, out);
	//Stormer_result[strlen(Stormer_result) - 3] = '\n';
	Stormer_result[strlen(Stormer_result) - 2] = '\0';
	printf("%s", Stormer_result);
	printf("\n");
	free(format_str);

	for(int i = 0; i < 4; i++) {
		mpfr_clear(arctg_f[i]);
		mpfr_clear(term_f[i]);
		mpfr_clear(frac_f[i]);
		mpz_clear(coef_z[i]);
		mpq_clear(frac_q[i]);
	}
	return 0;
}

char *generate_format_str(char *label, long int prec) {
	char *format_str = malloc((strlen(label) + 64) * sizeof(char));
	format_str[0] = '\0';
	strcat(format_str, label);
	strcat(format_str, "%.");
	char prec_str[32];
	sprintf(prec_str, "%ld", prec);
	strcat(format_str, prec_str);
	strcat(format_str, "Rf");
	return format_str;
}

long int compare_results(mpfr_t Takano_result, mpfr_t Stormer_result) {
	long int decimal_places = 0;
	char K_Takano_result[prec + 5];
	char FCM_Stormer_result[prec + 5];
	K_Takano_result[0] = '\0';
	FCM_Stormer_result[0] = '\0';
	mpfr_sprintf(K_Takano_result, generate_format_str("", prec - prec_buff), Takano_result);
	mpfr_sprintf(FCM_Stormer_result, generate_format_str("", prec - prec_buff), Stormer_result);

	for(int i = 0; i < strlen(K_Takano_result); i++) {
		if(K_Takano_result[i] == FCM_Stormer_result[i]) decimal_places++;
		else break;
	}

	return decimal_places - 2;
}