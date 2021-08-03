// Решение/работа с системой дифференциальных уравнений с n переменными и n неизвестными с учетом невязок
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
//*********************************************************
namespace StandartHelperLibrary.MathHelper
{
    /// <summary>
    /// Решение/работа с системой дифференциальных уравнений с n переменными и n неизвестными с учетом невязок
    /// </summary>
    public partial class TDifferentialSolver
    {
        //------------------------------------------------------------
        /// <summary>
        /// Решение системы дифференциальных уравнений с n переменными и n неизвестными  методом Рунге-Кутты 4ого порядка 
        /// </summary>
        /// <param name="Equation">Решаемая система уравнений и настройки решателя</param>
        /// <returns>Результат решения</returns>
        public static TSystemResultDifferential SolveSystemResidualFourRungeKutta(ISystemDifferentialEquation Equation)
        {
            double X = Equation.Min_X;                              // Крайняя левая точка диапазона "х" 
            double h = Equation.Step;                               // Шаг сетки "h" 
            int t = Equation.Rounding;                              // Округление до нужного знака, после запятой 
            int NumberOfIterations = Equation.CountIterations;      // Количество итераций
            int NumberOfEquations = Equation.CountEquations;        // кол-во урвнений
            List<double> InitArray = Equation.InitArray;            // Начальные значения Y для системы
            TSystemResultDifferential ResultSystemDifferential = new TSystemResultDifferential();

            //объявляются списки в которых будут храниться значения, которые были вычисленны. 
            List<double[]> Ys = new List<double[]>();

            // Массивы коэффициентов участвующих в методе Рунге-Кутты
            double[] Coefs1 = new double[NumberOfEquations];// число 1-ыx коэф. метода по числу уравнений
            double[] Coefs2 = new double[NumberOfEquations];// число 2-ыx коэф. метода по числу уравнений
            double[] Coefs3 = new double[NumberOfEquations];// число 3-иx коэф. метода по числу уравнений
            double[] Coefs4 = new double[NumberOfEquations];// число 4-ыx коэф. метода по числу уравнений

            double[] Y2 = new double[NumberOfEquations]; // число переменных для 2-го коэф. включая независимую
            double[] Y3 = new double[NumberOfEquations];// число переменных для 3-го коэф. включая независимую
            double[] Y4 = new double[NumberOfEquations];// число переменных для 4-го коэф. включая независимую

            //копируем начальные значения игреков в массив, который будет использоваться для вычислений
            double[] Y = new double[NumberOfEquations];//Массив всех Y
            for (int k = 0; k < Equation.InitArray.Count; k++)
            {
                Y[k] = InitArray[k];
            }

            //Закидываем первые(начальные) значения в поинт и далее в резалт
            TPointSystemDifferential PointSystemDifferentialInitial = new TPointSystemDifferential
            {
                Result = new double[NumberOfEquations],
                IndexIteration = 0,
                X = X
            };
            for (int i = 0; i < Y.Length; i++)
            {
                PointSystemDifferentialInitial.Result[i] = Y[i];
            }
            ResultSystemDifferential.SystemPoints.Add(PointSystemDifferentialInitial);

            //Массив текущих значений, по которым мы вычисляем невязку.(отслеживаемые значения, которые при достижения граничного, завершают интегрирование)
            List<double[]> Values = new List<double[]>();//лист массивов отслеживаемых значений
            var TrackedValues = CalculateValuesOfTrackedVariables(X, Y);//лист класса свойств граничных значений и вычисления невязки
            double[] Values_arr = new double[TrackedValues.Count()];//массив отслеживаемых значений на конкретной итерации
            for (int i = 0; i < TrackedValues.Count(); i++)//запись новых значений в массив
            {
                Values_arr[i] = TrackedValues[i].CurrentValue;
            }
            Values.Add(Values_arr);//добавление нового массива в лист массивов отслеживаемых значений
            
            //Цикл интегрирования. На каждой итерации совершается 1 "Шаг" интегрирования
            for (int i = 0; i < NumberOfIterations; i++)
            {
                double[] Ys_arr = new double[Y.Length];//создаем новый массив, который заполняем значениями игреков пердыдущей итерации
                for (int j = 0; j < Y.Length; j++)
                {
                    Ys_arr[j] = Y[j];
                }
                Ys.Add(Ys_arr);//добавление массива значений игреков в лист, который хранит значений игреков на всех итерациях. для возможности отката назад

                (double[] Y, double[] Coefs1, double[] Coefs2, double[] Coefs3, double[] Coefs4) Result;//объявляется вне цикла do while чтобы не была локальной и можно было перекидывать значения после выхода из цикла
                bool Incorrect = new bool();
                bool TimeToStop = new bool();
                bool StepChanged = new bool();
                TimeToStop = false;// Меняется на "правда", когда наши отслеживаенмые значения попадают в граничные с указанной точностью и как следствие мы заканчиваем интегрирование.
                do
                {
                    StepChanged = false;  // Изменяется на "правла", если шаг был изменен. Требуется для перезаписи старых значений У в массив Y[], чтобы расчитать новые значения с новым шагом.
                    Incorrect = true;     // Пока "правда" цикл do while продолжает работу
                    //Вычисляем новые значения Y и Coeffs для шага h
                    Result = CalculateValuesOf_Y(X, Y, h, Equation);

                    //просчет контрольных значений и невязок, так же хранит точность, название 
                    TrackedValues = CalculateValuesOfTrackedVariables(X + h, Y);

                    if (i < 1)//не проверяем выход по невязкам на первой(0-ой) итерации (т.к. нет еще пары значений текущей итерации и предыдущей, чтобы вычислить дельты) 
                        break;
                    
                    double[] Deltas = new double[TrackedValues.Count()];//создается массив для записи в него значений "дельт" - изменение отслеживаемого значения между текущей итерации и предыдущей
                    for (int j = 0; j < TrackedValues.Count(); j++)
                    {
                        Deltas[j] = TrackedValues[j].CurrentValue - Values[Values.Count() - 1][j];//вычисление дельты конкретного отслеживаемого значения
                    }

                    //проверка на необходимость уменьшения шага
                    for (int j = 0; j < TrackedValues.Count(); j++)
                    {
                        if (Deltas[j] > 0)//если да, то приближение к граничному значению производится снизу вверх
                        {
                            if (TrackedValues[j].CurrentValue > TrackedValues[j].BoundaryValue + TrackedValues[j].Accuracy)
                            {
                                h = h / 2d;
                                StepChanged = true;
                                break;
                            }
                        }
                        else if (Deltas[j] < 0) //если да, то приближение к граничному значению производится сверху вниз 
                        {
                            if (TrackedValues[j].CurrentValue < TrackedValues[j].BoundaryValue - TrackedValues[j].Accuracy)
                            {
                                h = h / 2d;
                                StepChanged = true;
                                break;
                            }
                        }
                    }

                    if (StepChanged)//если шаг изменен, то надо заного проверить не перескакиваем ли мы с уменьшенным шагом, а так же вычислить игреки заного. А для вычисления игреков заного надо перезаписать старые игреки в массив
                    {
                        for (int j = 0; j < Y.Length; j++)
                        {
                            Y[j] = Ys[Ys.Count() - 1][j];
                        }
                    }
                    else
                    {
                        Incorrect = false;
                        //проверка, не вошло ли уже значение в область граниченрого значения, а значит не пора ли уже завершить интегрирование
                        for (int j = 0; j < TrackedValues.Count(); j++)
                        {
                            if (Math.Abs(TrackedValues[j].Residual) <= TrackedValues[j].Accuracy)
                            {
                                TimeToStop = true;
                                break;
                            }
                        }
                    }
                }
                while (Incorrect);
                

                Values_arr = new double[TrackedValues.Count()];//создается новый массив для записи в него текущих значений отслеживаемых переменных
                for (int j = 0; j < TrackedValues.Count(); j++)
                {
                    Values_arr[j] = TrackedValues[j].CurrentValue;//заполнение массива текущими значениями(для данного шага) отслеживаемых переменных
                }
                Values.Add(Values_arr);//добавляется в лист массивов с отслеживаемыми переменными

                //прибавляем шаг к Х
                X += h;

                //создаем новый поинт и записываем результаты вычислений в него
                TPointSystemDifferential PointSystemDifferential = new TPointSystemDifferential
                {
                    Result = new double[NumberOfEquations],
                    IndexIteration = i + 1,
                    Coeffs = new List<double[]> { Result.Coefs1, Result.Coefs2, Result.Coefs3, Result.Coefs4 }
                };
                //Записываем новые значения в поинт, а далее поинт в резалт
                PointSystemDifferential.X = X;
                for (int j = 0; j < Y.Length; j++)
                {
                    PointSystemDifferential.Result[j] = Result.Y[j];
                }
                ResultSystemDifferential.SystemPoints.Add(PointSystemDifferential);

                if (TimeToStop)//если условие выполняется, то мы заканчиваем интегрирование и возвращем собранные в ходе вычислений значения
                    return ResultSystemDifferential;
            }
            // Вернуть результат
            return ResultSystemDifferential;
        }
 //-----------------------------------------------------------------------------------------------------------------------------------------------
        /// <summary>
        /// Метод вычисления новых значений Y с шагом в h
        /// </summary>
        /// <param name="X">Старое значение Х (для Х + h будут вычисляться Y)</param>
        /// <param name="Y">Старые значения Y (на основании которых будут найдены новые значения Y)</param>
        /// <param name="h">Значение "шага"</param>
        /// <param name="Equation">Решаемая система уравнений и настройки решателя</param>
        /// <returns></returns>
        private static (double[] Y, double[] Coefs1, double[] Coefs2, double[] Coefs3, double[] Coefs4) CalculateValuesOf_Y(double X, double[] Y, double h, ISystemDifferentialEquation Equation)
        {
            int NumberOfEquations = Equation.CountEquations;

            double Kx2_3 = X + h / 2;
            double Kx4 = X + h;

            var Coefs1 = Equation.ComputeEquation(X, Y);

            // Находим значения переменных для второго коэф. 
            var Y2 = CalcYVolumeForCoefs(Y, Coefs1, NumberOfEquations, h, true);
            var Coefs2 = Equation.ComputeEquation(Kx2_3, Y2);

            // Находим значения переменных для третьго коэф.
            var Y3 = CalcYVolumeForCoefs(Y, Coefs2, NumberOfEquations, h, true);
            var Coefs3 = Equation.ComputeEquation(Kx2_3, Y3);

            // Находим значения переменных для 4 коэф.
            var Y4 = CalcYVolumeForCoefs(Y, Coefs3, NumberOfEquations, h, false);
            var Coefs4 = Equation.ComputeEquation(Kx4, Y4);

            // Находим новые значения переменных включая независимую    
            for (int k = 0; k < NumberOfEquations; k++)
            {
                Y[k] += (Coefs1[k] + 2 * (Coefs2[k] + Coefs3[k]) + Coefs4[k]) * h / 6d;
            }
            return (Y, Coefs1, Coefs2, Coefs3, Coefs4);
        }
 //--------------------------------------------------------------------------------------------------------------------------------------------------------------
        /// <summary>
        /// Вычисляет значение Y для вычисление коэффициентов
        /// </summary>
        /// <param name="Y">Старые значения игреков на основании, которых вычисляются коэф., а затем на основании коэф. вычисляются новые значения игреков</param>
        /// <param name="Coefs">Массив соответствующих коэффициентов</param>
        /// <param name="NumberOfEquations">Число уравнений. Отправляется для выполнения цикла соответствующее кол-во раз</param>
        /// <param name="This_is_2nd_or_3rd_Y">правда если при вызове название данной булиновой переменной верно</param>
        /// <returns></returns>
        private static double[] CalcYVolumeForCoefs(double[] Y, double[] Coefs, int NumberOfEquations, double h, bool This_is_2nd_or_3rd_Y)
        {
            double[] Y_OUT = new double[NumberOfEquations];
            double K234 = new double(); //коэффициент, который только и отличается при расчете У2/У3/У4
            if (This_is_2nd_or_3rd_Y)   //определяется значение этого коэфа
                K234 = 2d;
            else
                K234 = 1d;
            for (int i = 0; i < NumberOfEquations; i++)
            {
                Y_OUT[i] = Y[i] + h * Coefs[i] / K234;       //выполняется прямое предназначение метода
            }
            return Y_OUT;
        }
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------
        /// <summary>
        /// Вывести отладочную информацию в консоль и в файл если задано имя
        /// </summary>
        /// <param name="Result">Результат решения системы дифф. уравнений</param>
        /// <param name="FileName">Имя файла</param>
        public static void Debug(TSystemResultDifferential Result, string FileName = "")
        {
            // В файл
            if (FileName.Length > 0) File.WriteAllText(FileName, Result.ToString());
            // В консоль
            Console.WriteLine(Result.ToString());
        }
 //--------------------------------------------------------------------------------------------------------------------------------------------------------------------
        /// <summary>
        /// Простой пример системы дифференциальных уравнений  и ее решения 
        /// <returns>Результат решения</returns>
        public static TSystemResultDifferential Example_dN_Residual()
        {
            // Создаем систему уравнений, которая должна решаться и задаем ее параметры 
            ISystemDifferentialEquation Equation = new TEquation_dN()
            {
                Equation = new AEquation_dN((X, Y) =>
                {
                    double[] FunArray = new double[1];//кол-во элементов в массиве должно быть = кол-во уравнений
                    //---------------------------------------------------------------
                    //задаются уравнения
                    FunArray[0] = (X * X - 2 * Y[0]);
                    //FunArray[1] = (1d);
                    //FunArray[2] = (X + Y[0] + 3 * Y[1] + Y[2] + Y[3] + Y[4]);
                    //FunArray[3] = (5 * X + 2 * Y[0] + 3 * Y[1] + Y[2] + Y[3] + Y[4]);
                    //FunArray[4] = (2 * X + 2 * Y[0] + 3 * Y[1] + Y[2] + Y[3] + Y[4]);
                    //---------------------------------------------------------------

                    return FunArray;   // интегрируемая система
                }),
                InitArray = new List<double> { 1/*, 0, 1, 1, 1, */},
                CountIterations = 1000,
                Min_X = 0,
                Rounding = 3,
                Step = 0.1,
                CountEquations = 1
            };
            // Решаем
            return SolveSystemResidualFourRungeKutta(Equation);
        }
//-------------------------------------------------------------------------------------------------------------------------------------
        /// <summary>
        /// Вычисляет значение отслеживаемых переменных 
        /// </summary>
        /// <param name="X">Новое значение X</param>
        /// <param name="Y">Новое значение У вычесленное по </param>
        /// <returns></returns>
        private static List<TBoundaryValue> CalculateValuesOfTrackedVariables(double X, double[] Y)
        {
            List<TBoundaryValue> BoundaryValues = new List<TBoundaryValue>();
            TBoundaryValue R1 = new TBoundaryValue()
            {
                Name = "", //название граничного значения(опционально и можно не указывать. создано чтобы не забывать какое граничное значение тут записано
                BoundaryValue = 0.33d, //граничное значение, по которому завершается интегрирование
                CurrentValue = Y[0], //вычисление текущего значения, которое мы сравниваем с граничным
                Accuracy = 0.01d //точность граничного значения
            };

            TBoundaryValue R2 = new TBoundaryValue()
            {
                Name = "",
                BoundaryValue = 1d,
                CurrentValue = X,
                Accuracy = 0.01d
            };

            //в лист записываются отслеживаемые переменные
            //Residuals.Add(R1);
            BoundaryValues.Add(R2);

            return BoundaryValues;
        }
    }
}