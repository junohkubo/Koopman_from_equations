CC = cc
SRC = $(OBJ:%.o=%.c)
OBJ = dual_computation.o 
TARGET = app_dual_computation

$(TARGET):$(OBJ)
	$(CC) -o $(TARGET) $(OBJ) -O3 -lm

clean:
	rm -f $(TARGET) $(OBJ)
