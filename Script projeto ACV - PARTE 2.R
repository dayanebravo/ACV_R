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

# Verificando se a matriz A é quadrada (mesmo número de linhas e colunas)
if (nrow(A) == ncol(A) & det(A) != 0) {
  print("A matriz é invertível.")
} else {
  stop("A matriz não é invertível. Revise a matriz A")
}

##################### Inventário M ################################
cat("##################### Inventário M ################################", "\n")

# Cálculo do inventário M=BA^(-1)k
inversa_A = solve(A)
s = inversa_A  %*% k
M = B %*% s
# Imprimindo o resultado do inventário
cat(M, "\n")

########################## Impacto h ##############################
cat("########################## Impacto h ##############################", "\n")

# Cálculo do impacto h=Q.M
h <- Q %*% M
# Imprimindo o resultado do impacto
cat(h, "\n")

################### Normalização h_til ############################
cat("################### Normalização h_til ############################", "\n")

# Cálculo da normalização h_til=h.h_ponto
h_til <- h / h_ponto
# Imprimindo o resultado do normalização
cat(h_til, "\n")

############## Ponderação w #######################################
cat("############## Ponderação w #######################################", "\n")

# Cálculo da ponderação w_calc=w*h_til
w_calc <- sum(w * h_til)
# Imprimindo o resultado da ponderação
cat(w_calc, "\n")

############# Sensibilidade de S em relação à A ###################
cat("############# Sensibilidade de S em relação à A ###################", "\n")

# Número de linhas e colunas da inversa da A e de S
num_linhas_invA <- nrow(inversa_A)
num_colunas_invA <- ncol(inversa_A)
num_linhas_s <- nrow(s)

negativa_invA = -inversa_A #negativa da invA

#criando vetores para guardar os elementos de s e invA
elementos_negativa <- numeric()
elementos_Sj <- numeric()
indice_neg <- 1
indice_sj <- 1

# A ordem de i, j e k é importante, pois vai organizar as posições dos vetores
# Eles serão multiplicados após essas interações, fora do laço
# Essa mesma ideia serve para os cálculos futuros
for (i in 1:num_colunas_invA) {
  for (j in 1:num_linhas_s) {
    for (k in 1:num_linhas_invA) {
      #guarda o elemento da posição k,i no vetor
      elementos_negativa[indice_neg] <- negativa_invA[k, i] 
      indice_neg <- indice_neg + 1 # troca para o próximo índice do vetor
      
      elementos_Sj[indice_sj] <- s[j] #guarda o elemento da posição j no vetor
      indice_sj <- indice_sj + 1 # troca para o próximo índice do vetor
      
    }
  }
}

ska = elementos_negativa * elementos_Sj # fórmula: [-(A^(-1))_{ki}*s_j]
# Imprimindo o resultado da sensibilidade de s em relação à A
cat(ska, "\n")

############### Sensibilidade de G em relação à A #################
cat("############### Sensibilidade de G em relação à A #################", "\n")

lambda = B %*% inversa_A #cálculo do lambda
negativa_lambda = -lambda # seu negativo

num_linhas_lamb <- nrow(lambda)
num_colunas_lamb <- ncol(lambda)

#criando vetores para guardar os elementos de s e neg_lambda
elementos_negativa_lamb <- numeric()
elementos_Sj <- numeric() #sobrescreve o anterior o sj anterior (faremos sempre)
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

gka = elementos_negativa_lamb * elementos_Sj  # fórmula: [-\lamb_{ki}*s_j]
# Imprimindo o resultado da sensibilidade de g
cat(gka, "\n")

############ Sensibilidade de G com relação à B ###################
cat("############ Sensibilidade de G com relação à B ###################", "\n")

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
      
      # quando i=j preenchemos o delta com 1, o restante é 0
      if (k == i) {  
        delta[indice_delta] <- 1
      }
      else {
        delta[indice_delta] <- 0
      }
      indice_delta <- indice_delta + 1
      
      elementos_Sj[indice_sj] <- s[j]
      indice_sj <- indice_sj + 1
    }
  }
}

gkb = delta * elementos_Sj # fórmula: [s_j*\delta_{ik}]
# Imprimindo o resultado da sensibilidade de g com delta
cat(gkb, "\n")

################## Sensibilidade de H com relação à A ##########################
cat("############# Sensibilidade de H com relação à A ##################", "\n")

num_linhas_q <- nrow(Q)
num_colunas_q <- ncol(Q)

elementos_Sj <- numeric()
indice_sj <- 1
qkl <- numeric()
indice_qkl <- 1
lambda_li <- numeric()
indice_lambda_li <- 1

for (i in 1:num_colunas_lamb) {
  for (j in 1:num_linhas_s) {
    for (k in 1:num_linhas_q) {
      for (l in 1:num_colunas_q) {
        
        elementos_Sj[indice_sj] <- -s[j]
        indice_sj <- indice_sj + 1
        
        qkl[indice_qkl] <- Q[k,l]
        indice_qkl <- indice_qkl + 1
        
        lambda_li[indice_lambda_li] <- lambda[l,i]
        indice_lambda_li <- indice_lambda_li + 1
        
      } 
    }
  }
}

qkl_lambda = qkl * lambda_li  # primeira parte da conta
produto_sj_qkl_lambda = elementos_Sj * qkl_lambda # segunda parte da conta

# TERCEIRA PARTE DA CONTA:
# Os 32 elementos serão divididos em 8 segmentos
# cada segmento contém 4 elementos que serão somados entre si.

# De forma genéria, a quantidade de segmentos é 
# o tamanho de qkl dividido pelo número de colunas de q
num_segmentos <- length(qkl_lambda)/num_colunas_q 

if (length(qkl_lambda) %% num_segmentos == 0){ #verifica o tamanho do vetor
  
  #cada segmento tem 4 elementos a serem somados
  #O tamanho do segmento (4), será o tamanho de qkl dividido pelo nº de segmentos
  tamanho_segmento <- length(qkl_lambda)/num_segmentos 
  
  # separando os segmentos por tamanho, para fazer a soma em seguida
  # o split corta o vetor de 32 elementos em 8 segmentos, com 4 elementos cada
  segmentos <- split(produto_sj_qkl_lambda, rep(1:num_segmentos, each = tamanho_segmento))
  
  #finalmente, faremos asa somas dos 4 elementos de cada segmento
  hk_aij <- numeric()
  for (i in 1:num_segmentos) {
    # somando os elementos de cada segmento:
    hk_aij[i] <- sum(segmentos[[i]])
  }
  # Imprimindo o resultado da sensibilidade de h com relação à A
  cat(hk_aij, "\n")  # fórmula: [-s_j*\sum_l(q_{kl}*\lam_{li})]
}

###################### Sensibilidade de H com relação à B ######################
cat("############# Sensibilidade de H com relação à B ##################", "\n")

elementos_Sj <- numeric()
qki <- numeric()
indice_sj <- 1
indice_qki <- 1

for (i in 1:num_colunas_q) {
  for (j in 1:num_linhas_s) {
    for (k in 1:num_linhas_q) {
      
      elementos_Sj[indice_sj] <- s[j]
      indice_sj <- indice_sj + 1
      
      qki[indice_qki] <- Q[k,i]
      indice_qki <- indice_qki + 1
      
    } 
  }
}

hk_bij = qki * elementos_Sj  # fórmula: [q_{ki}*s_j]

# Imprimindo o resultado da sensibilidade de h com relação à B
cat(hk_bij, "\n")
qki_sj = hk_bij  # eu troquei o nome, para depois padronizar as últimas contas

################## Sensibilidade de H com relação à Q ##########################
cat("############# Sensibilidade de H com relação à Q ##################", "\n")

num_linhas_g <- nrow(g)

gj <- numeric()
elementos_delta <- numeric()
indice_gj <- 1
indice_delta <- 1

for (i in 1:num_linhas_q) {
  for (j in 1:num_linhas_g) {
    for (k in 1:num_linhas_q) {
      
      # quando i=j preenchemos o delta com 1, o restante é 0
      if (k == i) {  
        elementos_delta[indice_delta] <- 1
      }
      else {
        elementos_delta[indice_delta] <- 0
      }
      indice_delta <- indice_delta + 1
      
      gj[indice_gj] <- g[j]
      indice_gj <- indice_gj + 1
      
    }
  }
}

hk_qij <- gj * elementos_delta  # fórmula: [g_j*\delta_{ik}]

# Imprimindo o resultado da sensibilidade de h com relação à Q
cat(hk_qij, "\n")

############################### INCERTEZA DE H #################################
cat('#################### INCERTEZA DE H ###############################', "\n")

# primeira parte da conta
#fórmula: [\sum_{i,j}(-s_j*\sum_l(q_{kl}*\lam_{li}))^2*var(a_{ij)]
num_linhas_h <- nrow(h)
indice_aux <- num_linhas_h - 1

elementos_varA <- numeric()
indice_varA <- 1

for (i in 1:num_colunas_lamb) {
  for (j in 1:num_linhas_s) {
    for (aux in 0:indice_aux) {
      
      elementos_varA[indice_varA+aux] <- varA[i,j]
    }
    indice_varA <- indice_varA + num_linhas_h
  }
}

calc1 <- ((hk_aij)**2) * elementos_varA 
var_h_aij <- calc1

# segunda parte da conta
# fórmula: [\sum_{i,j}(q_{ki}*s_j)^2*var(b_{ij)]

num_linhas_h <- nrow(h)
indice_aux <- num_linhas_h - 1

elementos_varB <- numeric()
indice_varB <- 1

for (i in 1:num_colunas_q) {
  for (j in 1:num_linhas_s) {
    for (aux in 0:indice_aux) {
      
      elementos_varB[indice_varB+aux] <- varB[i,j]
    }
    indice_varB <- indice_varB + num_linhas_h
  }
}

calc2 <- ((qki_sj)**2) * elementos_varB
var_h_bij <- calc2

# terceira parte da conta
# fórmula: [\sum_{j}(g_j)^2*var(q_{kj})]

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
var_h_qij <- calc3


cat("######### Incerteza de h em A #############", "\n")
cat(var_h_aij, "\n")
cat("######### Incerteza de h em B #############", "\n")
cat(var_h_bij, "\n")
cat("######### Incerteza de h em Q #############", "\n")
cat(var_h_qij, "\n")

######################## GRÁFICO DOS PONTOS CHAVES #############################
cat('#################### GRÁFICO DOS PONTOS CHAVES ####################', "\n")

# definindo o tamanho do eixo x
tamanho_hk_aij <- length(hk_aij)
tamanho_hk_bij <- length(hk_bij)
tamanho_hk_qij <- length(hk_qij)
tamanho_eixo_x <- tamanho_hk_aij + tamanho_hk_bij + tamanho_hk_qij

# criando uma matriz nula, com duas colunas, uma para cada eixo, 
# e numero de linhas = tamanho_eixo_x
eixos_grafico <- matrix(data = NA, nrow = tamanho_eixo_x, ncol=2)
n <- 1

# inserindo na matriz os valores da sensilibidade na primeira coluna
# e os valores da incerteza na segunda coluna
# para A, B e Q
for (x in 1:tamanho_hk_aij) {
  eixos_grafico[n,1] <- hk_aij[x]
  eixos_grafico[n,2] <- calc1[x]
  n <- n + 1
}
for (x in 1:tamanho_hk_bij) {
  eixos_grafico[n,1] <- hk_bij[x]
  eixos_grafico[n,2] <- calc2[x]
  n <- n + 1
}
for (x in 1:tamanho_hk_qij) {
  eixos_grafico[n,1] <- hk_qij[x]
  eixos_grafico[n,2] <- calc3[x]
  n <- n + 1
}

# dividindo a matriz inicial pela quantidade k do h
matrizes_resultantes <- list()
for (z in 1:num_linhas_h) {
  indices <- seq(z, nrow(eixos_grafico), by = num_linhas_h)
  matriz_resultante <- eixos_grafico[indices, ]
  matrizes_resultantes[[z]] <- matriz_resultante
}

# construindo o gráfico para cada uma das k matrizes de h
for (z in 1:num_linhas_h) {
  plot(matrizes_resultantes[[z]][, 1], 
       matrizes_resultantes[[z]][, 2],
       main = paste("Sensibilidade x Incerteza - h", z),
       xlab = "Sensibilidade", ylab = "Incerteza")

  # formatação do gráfico  
  text(matrizes_resultantes[[z]][, 1], 
       matrizes_resultantes[[z]][, 2],
       pos = 2,
       pch = 16,              # Tipo de símbolo (ponto preto)
       col = "blue"           # Cor dos pontos
       )
}


print('fim')

#labels = paste("Ponto", 1:nrow(matrizes_resultantes[[z]])), 



##### Inicializar variáveis para a soma de elementos em POSIÇÕES pares e ímpares #####
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

cat("######### Incerteza de h_1 #############", "\n")
cat(var_hk1, "\n")
cat("######### Incerteza de h_2 #############", "\n")
cat(var_hk2, "\n")




######################## GRÁFICO DOS PONTOS CHAVES #############################
cat('#################### GRÁFICO DOS PONTOS CHAVES ####################', "\n")

######### Construindo o eixo x para o gráfico ##################################
h1_eixo_x <- c()
h2_eixo_x <- c()

# Loop para montar o vetor de h_k, que será o eixo x - primeira parte (aij)
for (i in 1:length(hk_aij)) {
  if (i %% 2 == 0) {
    # Se o índice for par, adicione o elemento ao vetor, na última posição
    h2_eixo_x <- c(h2_eixo_x, hk_aij[i])
  } else {
    # Se o índice for ímpar, adicione o elemento ao vetor, na última posição
    h1_eixo_x <- c(h1_eixo_x, hk_aij[i])
  }
}

# Loop para montar o vetor de h_k, que será o eixo x - segunda parte (bij)
for (i in 1:length(hk_bij)) {
  if (i %% 2 == 0) {
    # Se o índice for par, adicione o elemento ao vetor, na última posição
    h2_eixo_x <- c(h2_eixo_x, hk_bij[i])
  } else {
    # Se o índice for ímpar, adicione o elemento ao vetor, na última posição
    h1_eixo_x <- c(h1_eixo_x, hk_bij[i])
  }
}

# Loop para montar o vetor de h_k, que será o eixo x - terceira parte (qij)
for (i in 1:length(hk_qij)) {
  if (i %% 2 == 0) {
    # Se o índice for par, adicione o elemento ao vetor, na última posição
    h2_eixo_x <- c(h2_eixo_x, hk_qij[i])
  } else {
    # Se o índice for ímpar, adicione o elemento ao vetor, na última posição
    h1_eixo_x <- c(h1_eixo_x, hk_qij[i])
  }
}

######### Construindo o eixo y para o gráfico ##################################
h1_eixo_y <- c()
h2_eixo_y <- c()

# Loop para montar o vetor de h_k, que será o eixo y - primeira parte (aij)
for (i in 1:length(calc1)) {
  if (i %% 2 == 0) {
    # Se o índice for par, adicione o elemento ao vetor, na última posição
    h2_eixo_y <- c(h2_eixo_y, calc1[i])
  } else {
    # Se o índice for ímpar, adicione o elemento ao vetor, na última posição
    h1_eixo_y <- c(h1_eixo_y, calc1[i])
  }
}

# Loop para montar o vetor de h_k, que será o eixo y - segunda parte (bij)
for (i in 1:length(calc2)) {
  if (i %% 2 == 0) {
    # Se o índice for par, adicione o elemento ao vetor, na última posição
    h2_eixo_y <- c(h2_eixo_y, calc2[i])
  } else {
    # Se o índice for ímpar, adicione o elemento ao vetor, na última posição
    h1_eixo_y <- c(h1_eixo_y, calc2[i])
  }
}

# precisamos inserir 4 valores nulos no eixo y para h2 antes do laço
valor_zero <- 0
h2_eixo_y <- append(h2_eixo_y, rep(valor_zero, num_colunas_q))

# Loop para montar o vetor de h_k, que será o eixo y - terceira parte (qij)
for (i in 1:length(calc3)) {
  if (i %% 2 == 0) {
    # Se o índice for par, adicione o elemento ao vetor, na última posição
    h2_eixo_y <- c(h2_eixo_y, calc3[i])
  } else {
    # Se o índice for ímpar, adicione o elemento ao vetor, na última posição
    h1_eixo_y <- c(h1_eixo_y, calc3[i])
  }
}

# precisamos inserir 4 valores nulos no eixo y para h1 depois do laço
h1_eixo_y <- append(h1_eixo_y, rep(valor_zero, num_colunas_q))


################ Criando um gráfico de dispersão para h1 ######################
plot(h1_eixo_x, h1_eixo_y, 
     xlab = "Sensibilidade",       # Rótulo do eixo X
     ylab = "Incerteza",       # Rótulo do eixo Y
     main = "Gráfico de Dispersão para h1", # Título do gráfico
     pch = 16,              # Tipo de símbolo (ponto preto)
     col = "blue"           # Cor dos pontos
)


################ Criando um gráfico de dispersão para h2 ######################
plot(h2_eixo_x, h2_eixo_y, 
     xlab = "Sensibilidade",       # Rótulo do eixo X
     ylab = "Incerteza",       # Rótulo do eixo Y
     main = "Gráfico de Dispersão para h2", # Título do gráfico
     pch = 16,              # Tipo de símbolo (ponto preto)
     col = "red"           # Cor dos pontos
)
