## Кодер и декодер Витерби
#### О проекте
Программа написана на языке с++ в среде Microsoft Visual Studio 2022 с использованием GNUplot для вывода графика. Код реализует кодер и декодер Витерби для кодов с параметром r = 1 / k, где k - количество полиномов используемых при кодировании/декодировании в двоичном симметричном канале с вероятностью ошибки на бит.

#### О написанном коде
Класс кодера:
* Конструктор - для кодирования сообщения нам необходимо знать только порождающие многочлены представленный в виде десятичного числа (x 3 + 1 -> 9, x 3 + x 2 + x + 1 -> 15 и т.д.);
* Функция кодера - считывает посимвольно поступающие биты и при помощи многочленов кодирует их в сверточный код.

Двоичный симметричный канал с вероятностью ошибки на бит:
* С заданной вероятностью изменяет бит на противоположный с 0 на 1, с 1 на 0.

Класс декодера Витерби
* Конструктор - для кодирования сообщения нам необходимо знать только порождающие многочлены представленный в виде десятичного числа (x 3 + 1 -> 9, x 3 + x 2 + x + 1 -> 15 и т.д.). В зависимости от k есть ограничение на количество сохраняемых путей - max_l_res;
* Функция декодера - считывает посимвольно поступающие биты сверточного кода и при помощи многочленов декодирует их в выходную последовательность.

#### График 
Вывод осуществляется при помощи GNUplot с использованием библиотеки gnuplot-iostream.h (https://github.com/dstahlke/gnuplot-iostream)

#### Написанный код в main
Код реализует вывод графика вероятности ошибки декодирования сверточного кода на бит при вероятности ошибки от 0 до 1 с шагом 0.01, при количестве тестов равных 100, входном сообщении "000000000000000000000000000000" и порождающими многочленами x 3 + x 2 + 1, x 3 + x 2 + x + 1, x 3 + x + 1.
