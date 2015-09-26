#include<stdio.h>
int main(int argc, char **argv){
	FILE *entrada;
	FILE *saida;
	float reader;
	int n = 0;
	float sum = 0.0;
	char str[] = {"m
	ediasCores.txt"};

	entrada = fopen(argv[1],"r");
	if(entrada==NULL){
		printf("Could not open file! Usage ./media \"YourFile.txt\"\n");
		return -1;
	}

	//fscanf(entrada,"%.3f",&reader);
	while((fscanf(entrada,"%f",&reader))!=EOF){
		sum+=reader;
		n++;
	}
	fclose(entrada);

	saida = fopen(str,"a");
	if(saida==NULL){
		printf("Erro ao abrir o arquivo: %s\n",str);
		return -1;
	}
	fprintf(saida,"%.3f\n",sum/n);
	fclose(saida);

	return 0;
}