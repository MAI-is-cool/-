// Решение/работа с дифференциальными уравнениями
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
    /// Решение/работа с дифференциальными уравнениями
    /// </summary>
    public partial class TDifferentialSolver
    {
//------------------------------------------------------------
        /// <summary>
        /// Решение дифференицального уравнения методом Рунге-Кутты 4ого порядка 
        /// </summary>
        /// <param name="Equation">Решаемое уравнение и настройки решателя</param>
        /// <returns>Результат решения</returns>
        public static TResultDifferential Solve_FourRungeKutta(IDifferentialEquation Equation)
        {
            // Результат
            TResultDifferential ResultDifferential = new TResultDifferential();
            // Рабочие переменные
            double x = Equation.Min_X;               // Крайняя левая точка диапазона "х" 
            double y = Equation.Min_Y;               // Значение "у" в крайней левой точке 
            double h = Equation.Step;                // Шаг сетки "h" 
            int t = Equation.Rounding;               // Округление до нужного знака, после запятой 
            int n = Equation.CountIterations;        // Количество итераций  
            // Интегрируем
            for (int j = 1; j <= n; j++)
            {
                TPointDifferential PointDifferential = new TPointDifferential();
                //
                double k1 = Equation.ComputeEquation(x, y);
                double k2 = Equation.ComputeEquation(x + h / 2, y + (h * k1 / 2));
                double k3 = Equation.ComputeEquation(x + h / 2, y + (h * k2 / 2));
                double k4 = Equation.ComputeEquation(x + h, y + (h * k3));
                //
                PointDifferential.Koeffs.Add("k1", Math.Round(k1, t));
                PointDifferential.Koeffs.Add("k2", Math.Round(k2, t));
                PointDifferential.Koeffs.Add("k3", Math.Round(k3, t));
                PointDifferential.Koeffs.Add("k4", Math.Round(k4, t));
                //
                double g = (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
                y = (g + y);
                //
                PointDifferential.Values.Add("x", Math.Round(x, t));
                PointDifferential.IndexIteration = j;
                PointDifferential.G = Math.Round(g, t);
                PointDifferential.Y = Math.Round(y, t);
                //
                ResultDifferential.Points.Add(PointDifferential);
                //
                x = x + h;
            }
            // Вернуть результат
            return ResultDifferential;
        }
//------------------------------------------------------------
        /// <summary>
        /// Вывести отладочную информацию в консоль и в файл если задано имя
        /// </summary>
        /// <param name="Result">Результат решения дифф. уравнения</param>
        /// <param name="FileName">Имя файла</param>
        public static void Debug(TResultDifferential Result, string FileName = "")
        {
            // В файл
            if (FileName.Length > 0) File.WriteAllText(FileName, Result.ToString());
            // В консоль
            Console.WriteLine(Result.ToString());
        }
//------------------------------------------------------------
        /// <summary>
        /// Простой пример дифференциального уравнения и его решения
        /// </summary>
        /// <returns>Результат решения</returns>
        public static TResultDifferential Example_dYdX()
        {
            // Создаем уравнение которое должно решаться и задаем ве параметры 
            IDifferentialEquation Equation = new TEquation_dYdX()
           
            {
                Equation = new AEquation_dYdX((X, Y) =>
                {
                    return X*X-Y;   // интегрируемая функция
                }),
                CountIterations = 10,
                Min_X = 0,
                Min_Z =1,
                Rounding = 3,
                Step = 0.1
            };
            
            // Решаем
            return Solve_FourRungeKutta(Equation);
        }
//------------------------------------------------------------
  
    }
}
