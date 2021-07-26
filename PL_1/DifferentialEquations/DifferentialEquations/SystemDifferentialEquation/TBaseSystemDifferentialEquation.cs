// Основа универсальной системы дифференциальных уравнений с n переменными и n неизвестными для последующего решения
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
//*********************************************************
namespace StandartHelperLibrary.MathHelper
{
    /// <summary>
    /// Основа универсальной системы дифференциальных уравнений с n переменными и n неизвестными для последующего решения
    /// </summary>
    public abstract class TBaseSystemDifferentialEquation
    {
        /// <summary>
        /// Крайняя левая точка диапазона "х" 
        /// </summary>
        public double Min_X { get; set; } = 0;

        /// <summary>
        /// Шаг интегрирования
        /// </summary>
        public double Step { get; set; } = 0.1;
        /// <summary>
        /// Количество итераций c шагом Step  
        /// </summary>
        public int CountIterations { get; set; } = 10;
        /// <summary>
        /// Округление до нужного знака, после запятой
        /// </summary>
        public int Rounding { get; set; } = 3;
        /// <summary>
        /// Количество уравнений системы 
        /// </summary>
        public int CountEquations { get; set; } = 5;

    }
}
