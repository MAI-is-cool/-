// Унифицированный интерфейс  системы дифференциальных уравнений с n переменными и n неизвестными 
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
//*********************************************************
namespace StandartHelperLibrary.MathHelper
{
    /// <summary>
    /// Унифицированный интерфейс  системы дифференциальных уравнений с n переменными и n неизвестными 
    /// </summary>
    public interface ISystemDifferentialEquation 
    {
        /// <summary>
        /// Начальные значения Y
        /// </summary>
        List<double> InitArray { get; set; }
        /// <summary>
        /// Крайняя левая точка диапазона "х" 
        /// </summary>
        double Min_X { get; set; }
        /// <summary>
        /// Шаг интегрирования
        /// </summary>
        double Step { get; set; }
        /// <summary>
        /// Количество итераций c шагом Step  
        /// </summary>
        int CountIterations { get; set; }
        /// <summary>
        /// Округление до нужного знака, после запятой
        /// </summary>
        int Rounding { get; set; }
        /// <summary>
        /// Рассчитать хранимое уравнение
        /// </summary>
        /// <param name="X">Переменная, по которой идет интегриривание</param>
        /// <param name="Y">Массив неизветных функций, зависящих от x </param>
        /// <returns>Результат</returns>
        double[] ComputeEquation(double X, double [] Y);
        /// <summary>
        /// Количество уравнений системы
        /// </summary>
        int CountEquations { get; set; }
    }
}

