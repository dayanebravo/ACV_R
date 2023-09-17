
# install.packages("readxl")
# Carregar a biblioteca readxl
library(readxl)

################################################

# Ler os arquivos XLSX com as matrizes
caminho_A <- "A.xlsx"
dados_A <- read_xlsx(caminho_A)
# remover a primeira coluna
dados_A <- dados_A[-1]
# converter df em matriz
A <- as.matrix(dados_A)

# fazendo o mesmo para as outras entradas
caminho_B <- "B.xlsx"
dados_B <- read_xlsx(caminho_B)
dados_B <- dados_B[-1]
B <- as.matrix(dados_B)

caminho_Q <- "Q.xlsx"
dados_Q <- read_xlsx(caminho_Q)
dados_Q <- dados_Q[-1]
Q <- as.matrix(dados_Q) 

caminho_h <- "hponto.xlsx"
dados_h <- read_xlsx(caminho_h)
dados_h <- dados_h[-1]
h_ponto <- as.matrix(dados_h)

caminho_w <- "w.xlsx"
dados_w <- read_xlsx(caminho_w)
dados_w <- dados_w[-1]
w <- as.matrix(dados_w)

caminho_k <- "f.xlsx"
dados_k <- read_xlsx(caminho_k)
dados_k <- dados_k[-1]
k <- as.matrix(dados_k)
f <- as.matrix(dados_k)


caminho_varA <- "varA.xlsx"
dados_varA <- read_xlsx(caminho_varA)
dados_varA <- dados_varA[-1]
varA <- as.matrix(dados_varA)

caminho_varB <- "varB.xlsx"
dados_varB <- read_xlsx(caminho_varB)
dados_varB <- dados_varB[-1]
varB <- as.matrix(dados_varB)

caminho_varQ <- "varQ.xlsx"
dados_varQ <- read_xlsx(caminho_varQ)
dados_varQ <- dados_varQ[-1]
varQ <- as.matrix(dados_varQ)

################################################


# Verificando se a matriz é quadrada (mesmo número de linhas e colunas)

if (nrow(A) == ncol(A) & det(A) != 0) {
  print("A matriz é invertível.")
} else {
  stop("A matriz não é invertível. Revise a matriz A")
}

################################################

# Cálculo do inventário M=BA^(-1)k
inversa_A = solve(A)

p1 <- B %*% inversa_A
M = p1  %*% k

# Imprimindo o resultado do inventário
print(M)

################################################

# Cálculo do impacto h=Q.M
h <- Q %*% M

# Imprimindo o resultado do impacto
print(h)

#############################################

# Cálculo da normalização h_til=h.h_ponto
h_til <- h / h_ponto

# Imprimindo o resultado do normalização
print(h_til)


################################################


# Cálculo da ponderação w_calc=w*h_til
w_calc <- sum(w * h_til)

# Imprimindo o resultado do ponderação
print(w_calc)