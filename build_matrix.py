from random import randint, uniform


size = int(input("Tamanho do sistema:"))

f = open("matrix"+ str(size) + ".mtx", "w")
f.write(str(size)+"\n")

for i in range(size):
    p = size * 2
    main = randint(1, p)
    sum = main
    for j in range(size):
        if i == j:
            string = " ".join((str(i + 1), str(j + 1), str(main), '\n'))
            f.write(string)            
        else:
            chute = uniform(0, 1)
            if chute >= .75:
                continue
            p = int(size / 2)
            value = randint(1, p)
            if sum - value > 0:
                sum = sum -value
                string = " ".join((str(i + 1), str(j + 1), str(value), '\n'))
                f.write(string)              
    value = uniform(0, 100)
    string = " ".join((str(i + 1), str(size + 1), str(value), '\n'))
    f.write(string)  
f.close()