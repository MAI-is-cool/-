// Решение/работа с дифференциальными уравнениями второго порядка 
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
        /// Решение линейного неоднородного диф. уравнения второго порядка
        /// </summary>
        /// <param name="Equation">Решаемое уравнение и настройки решателя</param>
        /// <returns>Результат решения</returns>
        public static TResultDifferentialTwo Solve_LinearInhomogeneousSecondOrderEquation(IDifferentialEquation Equation, IDifferentialEquation Equations)
        {
            // Результат
            TResultDifferentialTwo ResultDifferential = new TResultDifferentialTwo();
            // Рабочие переменные
            double x = Equation.Min_X;               // Крайняя левая точка диапазона "х" 
            double y = Equation.Min_Y;               // Значение "у" в крайней левой точке 
            double h = Equation.Step;                // Шаг сетки "h" 
            int t = Equation.Rounding;               // Округление до нужного знака, после запятой 
            int n = Equation.CountIterations;        // Количество итераций  
            double z = Equation.Min_Z;               // Значение "z" в крайней левой точке 
            for (int j = 1; j <= n; j++)
            {
                TPointDifferentialTwo PointDifferential = new TPointDifferentialTwo();
                //
                double k1 = Equation.ComputeEquation(x, y, z);
                double m1 = Equations.ComputeEquation(x, y, z);
                double k2 = Equation.ComputeEquation(x + h / 2, y + (k1 / 2), z + (m1 / 2));
                double m2 = Equations.ComputeEquation(x + h / 2, y + (k1 / 2), z + (m1 / 2));
                double k3 = Equation.ComputeEquation(x + h / 2, y + (k2 / 2), z + (m2 / 2));
                double m3 = Equations.ComputeEquation(x + h / 2, y + (k2 / 2), z + (m2 / 2));
                double k4 = Equation.ComputeEquation(x + h, y + k3, z + m3);
                double m4 = Equations.ComputeEquation(x + h, y + k3, z + m3);
                PointDifferential.Koeffs.Add("k1", Math.Round(k1, t));
                PointDifferential.Koeffs.Add("k2", Math.Round(k2, t));
                PointDifferential.Koeffs.Add("k3", Math.Round(k3, t));
                PointDifferential.Koeffs.Add("k4", Math.Round(k4, t));
                PointDifferential.Koeffs.Add("m1", Math.Round(m1, t));
                PointDifferential.Koeffs.Add("m2", Math.Round(m2, t));
                PointDifferential.Koeffs.Add("m3", Math.Round(m3, t));
                PointDifferential.Koeffs.Add("m4", Math.Round(m4, t));
                //
                double g = (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
                double l = (h / 6) * (m1 + 2 * m2 + 2 * m3 + m4);
                z = (l + z);
                y = (g + y);
               
                PointDifferential.Values.Add("x", Math.Round(x, t));
                PointDifferential.IndexIteration = j;
                PointDifferential.G = Math.Round(g, t);
                PointDifferential.L = Math.Round(l, t);
                PointDifferential.Z = Math.Round(z, t);
                PointDifferential.Y = Math.Round(y, t);
                
                ResultDifferential.PointsTwo.Add(PointDifferential);
                //
                x = x + h;

            }
            return ResultDifferential;
        }
//------------------------------------------------------------
        /// <summary>
        /// Вывести отладочную информацию в консоль и в файл если задано имя
        /// </summary>
        /// <param name="Result">Результат решения дифф. уравнения</param>
        /// <param name="FileName">Имя файла</param>
        public static void Debug(TResultDifferentialTwo Result, string FileName = "")
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
        public static TResultDifferentialTwo Example_dXdYdZ()
        {
            // Создаем уравнение которое должно решаться и задаем ве параметры 
            IDifferentialEquation Equation = new TEquation_dXdYdZ()

            {
                Equation = new AEquation_dXdYdZ((X, Y, Z) =>
                {
                    return X+ Y + Z;   // интегрируемая функция
                }),
                CountIterations = 10,
                Min_X = 0,
                Min_Y = 1,
                Min_Z = 1,
                Rounding = 3,
                Step = 1
            };
            IDifferentialEquation Equations = new TEquation_dXdYdZ()
            {
                Equation = new AEquation_dXdYdZ((X, Y, Z) =>
                {
                    return X + 2*Y+ Z;
                }),
                CountIterations = 10,
                Min_X = 0,
                Min_Y = 1,
                Min_Z = 1,
                Rounding = 3,
                Step = 1
            };
            // Решаем
            return Solve_LinearInhomogeneousSecondOrderEquation(Equation, Equations);
        }
        //------------------------------------------------------------
    }
}