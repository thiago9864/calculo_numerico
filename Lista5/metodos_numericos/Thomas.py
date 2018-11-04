class Thomas():

    def executar(self, M, d):
        ''' Resolve Ax = d onde A e uma matriz tridiagonal composta pelos vetores a, b, c
            a - subdiagonal
            b - diagonal principal
            c - superdiagonal.
        Retorna x
        '''
        passos = 0
        num_particoes = len(M[0])
        a = [] #subdiagonal
        b = [] #diagonal principal
        c = [] #superdiagonal
        
        #separa os vetores da matriz
        for i in range(num_particoes):
            for j in range(num_particoes):
                passos += 1
                if(i==j):
                    #preenche diagonal principal
                    b.append(M[i][j])
                elif(i==j+1):
                    #preenche diagonal inferior
                    a.append(M[i][j])
                elif(i==j-1):
                    #preenche diagonal superior
                    c.append(M[i][j])


        n = len(d) # len(d) == len(b)
        c_ = [ c[0] / b[0] ]
        d_ = [ d[0] / b[0] ]
        
        for i in range(1, n):
            passos += 1
            aux = b[i] - c_[i-1]*a[i-1]
            if i < n-1:
                c_.append( c[i] / aux )
            d_.append( (d[i] - d_[i-1]*a[i-1])/aux )
        
        # Substituicao de volta
        x = [d_[-1]]
        for i in range(n-2, -1, -1):
            passos += 1
            x = [ d_[i] - c_[i]*x[0] ] + x
        
        return [x, passos]