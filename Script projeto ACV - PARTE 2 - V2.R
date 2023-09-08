# caso a biblioteca não exista no pc
# install.packages("readxl")
# Carregar a biblioteca readxl
library(readxl)

################### LEITURA DAS ENTRADAS EM XLSX ###############################
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

###############################################################################

# Verificando se a matriz é quadrada (mesmo número de linhas e colunas)
if (nrow(A) == ncol(A) & det(A) != 0) {
  print("A matriz é invertível.")
} else {
  stop("A matriz não é invertível. Revise a matriz A")
}

###############################################################################

# Cálculo do inventário M=BA^(-1)k
inversa_A = solve(A)
s = inversa_A  %*% k
M = B %*% s
# Imprimindo o resultado do inventário
print(M)

###############################################################################

# Cálculo do impacto h=Q.M
h <- Q %*% M
# Imprimindo o resultado do impacto
print(h)

###############################################################################

# Cálculo da normalização h_til=h.h_ponto
h_til <- h / h_ponto
# Imprimindo o resultado do normalização
print(h_til)

###############################################################################

# Cálculo da ponderação w_calc=w*h_til
w_calc <- sum(w * h_til)
# Imprimindo o resultado da ponderação
print(w_calc)

########################### SENSIBILIDADE DE S ################################

# Número de linhas e colunas da inversa da A e de S
num_linhas_invA <- nrow(inversa_A)
num_colunas_invA <- ncol(inversa_A)
num_linhas_s <- nrow(s)

negativa_invA = -inversa_A #negativa da invA

#criando vetores para guardar os elementros de s e invA
elementos_negativa <- numeric()
elementos_Sj <- numeric()
indice_neg <- 1
indice_sj <- 1

#a ordem de i, j e k é importante, pois vai organizar as posições dos vetores
#eles serão multiplicados após essas interações
for (i in 1:num_colunas_invA) {
  for (j in 1:num_linhas_s) {
    for (k in 1:num_linhas_invA) {
      elementos_negativa[indice_neg] <- negativa_invA[k, i] #guarda o elemento da posição k,i no vetor
      indice_neg <- indice_neg + 1
      elementos_Sj[indice_sj] <- s[j] #guarda o elemento da posição j no vetor
      indice_sj <- indice_sj + 1
    }
  }
}

ska = elementos_negativa * elementos_Sj
# Imprimindo o resultado da sensibilidade de s
print(ska)
print('######################################################################')
###################### SENSIBILIDADE DE G #####################################

lambda = B %*% inversa_A #cálculo do lambda

num_linhas_lamb <- nrow(lambda)
num_colunas_lamb <- ncol(lambda)

negativa_lambda = -lambda

#criando vetores para guardar os elementros de s e neg_lambda
elementos_negativa_lamb <- numeric()
elementos_Sj <- numeric() #novo vetor para sj, sobrescreve o anterior
indice_neg <- 1
indice_sj <- 1

for (i in 1:num_colunas_lamb) {
  for (j in 1:num_linhas_s) {
    for (k in 1:num_linhas_lamb) {
      elementos_negativa_lamb[indice_neg] <- negativa_lambda[k, i]
      indice_neg <- indice_neg + 1
      elementos_Sj[indice_sj] <- s[j]
      indice_sj <- indice_sj + 1
    }
  }
}

gka = elementos_negativa_lamb * elementos_Sj
# Imprimindo o resultado da sensibilidade de g
print(gka)
print('######################################################################')
###################### SENSIBILIDADE DE G COM DELTA ##########################

g = M #já havíamos calculado, porém com outra notação M = B %*% s
num_linhas_b <- nrow(B)
num_linhas_g <- nrow(g)

elementos_Sj <- numeric()
delta <- numeric()
indice_delta <- 1
indice_sj <- 1

for (i in 1:num_linhas_b) {
  for (j in 1:num_linhas_s) {
    for (k in 1:num_linhas_g) {
      if (k == i) {  # quando i=j preenchemos o delta com 1, o restante é 0
        delta[indice_delta] <- 1
      }
      else {
        delta[indice_delta] <- 0
      }

      elementos_Sj[indice_sj] <- s[j]
      indice_sj <- indice_sj + 1
      indice_delta <- indice_delta + 1
    }
  }
}

gkb = delta * elementos_Sj
# Imprimindo o resultado da sensibilidade de g com delta
print(gkb)
print('######################################################################')
###################### SENSIBILIDADE DE H parte 1 ##########################

num_linhas_q <- nrow(Q)
num_colunas_q <- ncol(Q)

elementos_Sj <- numeric()
elementos_varA <- numeric()
qkl <- numeric()
lambda_li <- numeric()
indice_lambda_li <- 1
indice_sj <- 1
indice_qkl <- 1
indice_varA <- 1

for (i in 1:num_colunas_lamb) {
  for (j in 1:num_linhas_s) {
    for (k in 1:num_linhas_q) {
      for (l in 1:num_colunas_q) {
        elementos_Sj[indice_sj] <- -s[j]
        qkl[indice_qkl] <- Q[k,l]
        lambda_li[indice_lambda_li] <- lambda[l,i]
        
        indice_sj <- indice_sj + 1
        indice_lambda_li <- indice_lambda_li + 1
        indice_qkl <- indice_qkl + 1
      } 
    }
    # organizando os dados da variância para futuramente calcular a incerteza
    elementos_varA[indice_varA] <- varA[i,j]
    elementos_varA[indice_varA + 1] <- varA[i,j]
    indice_varA <- indice_varA + 2
  }
}

qkl_lambda = qkl * lambda_li
print(qkl_lambda)
print('######################################################################')
###################### SENSIBILIDADE DE H parte 2 ##########################

produto_sj_qkl_lambda = elementos_Sj*qkl_lambda

num_segmentos <- length(qkl_lambda)/num_colunas_q #os 32 elementos serão divididos em 8 segmentos

if (length(qkl_lambda) %% num_segmentos==0){ #verifica se o tamanho do vetor está certo
  tamanho_segmento <- length(qkl_lambda)/num_segmentos #cada segmento tem 4 elementos a serem somados
  
  # separando os segmentos por tamanho, para fazer a soma em seguida:
  segmentos <- split(produto_sj_qkl_lambda, rep(1:num_segmentos, each = tamanho_segmento))
  
  hk_aij <- numeric()
  for (i in 1:num_segmentos) {
    # somando os elementos de cada segmento:
    hk_aij[i] <- sum(segmentos[[i]])
  }
  # Imprimindo o resultado da sensibilidade de h
  print(hk_aij)
}
print('######################################################################')
###################### SENSIBILIDADE DE H parte 3 ##########################

elementos_Sj <- numeric()
elementos_varB <- numeric()
qki <- numeric()
indice_sj <- 1
indice_qki <- 1
indice_varB <- 1

for (i in 1:num_colunas_q) {
  for (j in 1:num_linhas_s) {
    for (k in 1:num_linhas_q) {
      elementos_Sj[indice_sj] <- s[j]
      qki[indice_qki] <- Q[k,i]

      indice_sj <- indice_sj + 1
      indice_qki <- indice_qki + 1
    } 
    # organizando os dados da variância para futuramente calcular a incerteza
    elementos_varB[indice_varB] <- varB[i,j]
    elementos_varB[indice_varB + 1] <- varB[i,j]
    indice_varB <- indice_varB + 2
  }
}

hk_bij = qki * elementos_Sj
print(hk_bij)
print('######################################################################')
###################### SENSIBILIDADE DE H parte 4 ##########################

num_linhas_g <- nrow(g)

elementos_varQ <- numeric()
gj <- numeric()
elementos_delta <- numeric()
indice_gj <- 1
indice_delta <- 1
indice_varQ <- 1

for (i in 1:num_linhas_q) {
  for (j in 1:num_linhas_g) {
    for (k in 1:num_linhas_q) {
      if (k == i) {  # quando i=j preenchemos o delta com 1, o restante é 0
        elementos_delta[indice_delta] <- 1
      }
      else {
        elementos_delta[indice_delta] <- 0
      }
      gj[indice_gj] <- g[j]
      indice_gj <- indice_gj + 1
      indice_delta <- indice_delta + 1
    
      elementos_varQ[indice_varQ] <- varQ[k,j]
      indice_varQ <- indice_varQ + 1
    }
  }
}

hk_qij <- gj * elementos_delta
print(hk_qij)
print('######################################################################')
############################### INCERTEZA DE H #################################
print('#################### INCERTEZA DE H ###################################')

calc1 <- ((hk_aij)**2) * elementos_varA
calc2 <- ((qki_sj)**2) * elementos_varB

elementos_varQ <- numeric()
gj <- numeric()
indice_gj <- 1
indice_varQ <- 1

for (j in 1:num_linhas_g) {
  for (k in 1:num_linhas_q) {
    gj[indice_gj] <- g[j]
    indice_gj <- indice_gj + 1
    elementos_varQ[indice_varQ] <- varQ[k,j]
    indice_varQ <- indice_varQ + 1
  }
}


calc3 <- (gj**2)*elementos_varQ

# Inicializar variáveis para a soma de elementos em posições pares e ímpares
soma_pares_calc1 <- 0
soma_impares_calc1 <- 0
soma_pares_calc2 <- 0
soma_impares_calc2 <- 0
soma_pares_calc3 <- 0
soma_impares_calc3 <- 0

# Loop através do vetor para calcular as somas
for (i in 1:length(calc1)) {
  if (i %% 2 == 0) {
    # Se o índice for par, adicione o elemento à soma de pares
    soma_pares_calc1 <- soma_pares_calc1 + calc1[i]
  } else {
    # Se o índice for ímpar, adicione o elemento à soma de ímpares
    soma_impares_calc1 <- soma_impares_calc1 + calc1[i]
  }
}

# Loop através do vetor para calcular as somas
for (i in 1:length(calc2)) {
  if (i %% 2 == 0) {
    # Se o índice for par, adicione o elemento à soma de pares
    soma_pares_calc2 <- soma_pares_calc2 + calc2[i]
  } else {
    # Se o índice for ímpar, adicione o elemento à soma de ímpares
    soma_impares_calc2 <- soma_impares_calc2 + calc2[i]
  }
}

# Loop através do vetor para calcular as somas
for (i in 1:length(calc3)) {
  if (i %% 2 == 0) {
    # Se o índice for par, adicione o elemento à soma de pares
    soma_pares_calc3 <- soma_pares_calc3 + calc3[i]
  } else {
    # Se o índice for ímpar, adicione o elemento à soma de ímpares
    soma_impares_calc3 <- soma_impares_calc3 + calc3[i]
  }
}

var_hk1 <- soma_impares_calc1 + soma_impares_calc2 + soma_impares_calc3
var_hk2 <- soma_pares_calc1 + soma_pares_calc2 + soma_pares_calc3

print(var_hk1)
print(var_hk2)
######################## GRÁFICO DOS PONTOS CHAVES #############################
print('#################### GRÁFICO DOS PONTOS CHAVES ########################')

######### construindo o eixo x para o gráfico ##################################
h1_eixo_x <- c()
h2_eixo_x <- c()

# Loop através do vetor para calcular as somas
for (i in 1:length(hk_aij)) {
  if (i %% 2 == 0) {
    # Se o índice for par, adicione o elemento ao vetor
    h2_eixo_x <- c(h2_eixo_x, hk_aij[i])
  } else {
    # Se o índice for ímpar, adicione o elemento ao vetor
    h1_eixo_x <- c(h1_eixo_x, hk_aij[i])
  }
}

# Loop através do vetor para calcular as somas
for (i in 1:length(hk_bij)) {
  if (i %% 2 == 0) {
    # Se o índice for par, adicione o elemento ao vetor
    h2_eixo_x <- c(h2_eixo_x, hk_bij[i])
  } else {
    # Se o índice for ímpar, adicione o elemento ao vetor
    h1_eixo_x <- c(h1_eixo_x, hk_bij[i])
  }
}

# Loop através do vetor para calcular as somas
for (i in 1:length(hk_qij)) {
  if (i %% 2 == 0) {
    # Se o índice for par, adicione o elemento ao vetor
    h2_eixo_x <- c(h2_eixo_x, hk_qij[i])
  } else {
    # Se o índice for ímpar, adicione o elemento ao vetor
    h1_eixo_x <- c(h1_eixo_x, hk_qij[i])
  }
}


######### construindo o eixo y para o gráfico ##################################
h1_eixo_y <- c()
h2_eixo_y <- c()

# Loop através do vetor para calcular as somas
for (i in 1:length(calc1)) {
  if (i %% 2 == 0) {
    # Se o índice for par, adicione o elemento ao vetor
    h2_eixo_y <- c(h2_eixo_y, calc1[i])
  } else {
    # Se o índice for ímpar, adicione o elemento ao vetor
    h1_eixo_y <- c(h1_eixo_y, calc1[i])
  }
}

# Loop através do vetor para calcular as somas
for (i in 1:length(calc2)) {
  if (i %% 2 == 0) {
    # Se o índice for par, adicione o elemento ao vetor
    h2_eixo_y <- c(h2_eixo_y, calc2[i])
  } else {
    # Se o índice for ímpar, adicione o elemento ao vetor
    h1_eixo_y <- c(h1_eixo_y, calc2[i])
  }
}


#############################################################
############ colocar os zeros ##########################
###########################################################3
# Loop através do vetor para calcular as somas
valor_zero <- 0
h2_eixo_y <- append(h2_eixo_y, rep(valor_zero, num_colunas_q))

for (i in 1:length(calc3)) {
  if (i %% 2 == 0) {
    # Se o índice for par, adicione o elemento ao vetor
    h2_eixo_y <- c(h2_eixo_y, calc3[i])
  } else {
    # Se o índice for ímpar, adicione o elemento ao vetor
    h1_eixo_y <- c(h1_eixo_y, calc3[i])
  }
}

h1_eixo_y <- append(h1_eixo_y, rep(valor_zero, num_colunas_q))


# Criar um gráfico de dispersão para h1
plot(h1_eixo_x, h1_eixo_y, 
     xlab = "Sensibilidade",       # Rótulo do eixo X
     ylab = "Incerteza",       # Rótulo do eixo Y
     main = "Gráfico de Dispersão para h1", # Título do gráfico
     pch = 16,              # Tipo de símbolo (ponto preto)
     col = "blue"           # Cor dos pontos
)


# Criar um gráfico de dispersão para h2
plot(h2_eixo_x, h2_eixo_y, 
     xlab = "Sensibilidade",       # Rótulo do eixo X
     ylab = "Incerteza",       # Rótulo do eixo Y
     main = "Gráfico de Dispersão para h2", # Título do gráfico
     pch = 16,              # Tipo de símbolo (ponto preto)
     col = "red"           # Cor dos pontos
)