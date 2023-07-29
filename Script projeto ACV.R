# install.packages("readxl")
# Carregar a biblioteca readxl
library(readxl)

################################################


# Ler o arquivo XLSX
caminho_A <- "A.xlsx"
dados_A <- read_xlsx(caminho_A)
A <- as.matrix(dados_A)

caminho_B <- "B.xlsx"
dados_B <- read_xlsx(caminho_B)
B <- as.matrix(dados_B)

caminho_Q <- "Q.xlsx"
dados_Q <- read_xlsx(caminho_Q)
Q <- as.matrix(dados_Q)

caminho_h <- "hponto.xlsx"
dados_h <- read_xlsx(caminho_h)
h_ponto <- as.matrix(dados_h)

caminho_w <- "w.xlsx"
dados_w <- read_xlsx(caminho_w)
w <- as.matrix(dados_w)

caminho_k <- "k.xlsx"
dados_k <- read_xlsx(caminho_k)
k <- as.matrix(dados_k)

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

################################################


# Cálculo da normalização h_til=h.h_ponto
h_til <- h / h_ponto

# Imprimindo o resultado do normalização
print(h_til)


################################################


# Cálculo da ponderação w_calc=w*h_til
w_calc <- sum(w * h_til)

# Imprimindo o resultado do ponderação
print(w_calc)