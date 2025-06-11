/// @file DeconvolveTimerUtils.h
///
/// @copyright (c) 2022 CSIRO
/// Australia Telescope National Facility (ATNF)
/// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
/// PO Box 76, Epping NSW 1710, Australia
/// atnf-enquiries@csiro.au
///
/// This file is part of the ASKAP software distribution.
///
/// The ASKAP software distribution is free software: you can redistribute it
/// and/or modify it under the terms of the GNU General Public License as
/// published by the Free Software Foundation; either version 2 of the License,
/// or (at your option) any later version.
///
/// This program is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with this program; if not, write to the Free Software
/// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
///
/// @author Minh Vuong <minh.vvuong@csiro.au>
///

#ifndef ASKAP_DECONVOLVE_TIMER_UTILS_H
#define ASKAP_DECONVOLVE_TIMER_UTILS_H

#include <map>
#include <string>
#include <memory>
#include <ctime>
#include <mutex>

namespace askap {
namespace utils {

enum class State {START, STOP};

/// @brief An abstract class which defines the interface for the timer.
class ITimer {
    public:
        virtual ~ITimer() {}
        /// @brief - start the time
        virtual void start() = 0;
        /// @brief - stop the time
        virtual void stop() = 0;
        /// @return a string summary of the timer
        virtual std::string summary() const = 0;
        /// @return the timer elapsed time when the timer
        ///  first started and the last timer it is stopped.
        virtual double elapsedTime() const = 0;

    protected:
        /// brief - timer's internal state
        askap::utils::State itsState;
        std::mutex itsMutex;
};

/// @detail A general timer class which is accessable from client code.
///         It determines the real implementation timer to use based
///         on what is available on the system. The preferred timer it
///         uses is in this order (1) MPI_Wtime, (2) omp_get_wtime and
///         std::chrono::system_clock
class Timer final : public ITimer {
    public:
        Timer();
        Timer(const Timer&) = delete;
        Timer& operator=(const Timer&) = delete; 
        ~Timer() {}
        /// @brief delegate the call to itsTimerImpl start
        void start() override;
        /// @brief delegate the call to itsTimerImpl stop
        void stop() override;
        /// @brief delegate the call to itsTimerImpl summary
        std::string summary() const override;
        /// @brief delegate the call to itsTimerImpl elapsedTime
        double elapsedTime() const override;

    private:
        /// @brief Real implementation timer.
        std::shared_ptr<ITimer> itsTimerImpl;
};

/// @brief A Standard timer class using the std::chrono::system_clock
class StdTimer final : public ITimer
{
  public:
    friend class Timer;
    ~StdTimer();
  private:
    StdTimer();
    StdTimer(const StdTimer&) = delete;
    StdTimer& operator=(const StdTimer&) = delete;

    /// @brief implementation of ITimer::start
    void start() override;
    /// @brief implementation of ITimer::stop
    void stop() override;
    /// @brief implementation of ITimer::summary
    double elapsedTime() const override;
    /// @brief implementation of ITimer::elapsedTime
    std::string summary() const override;

    std::time_t itsElapsedTime;
    std::time_t itsStartTime;
    std::time_t itsStopTime;
};

#ifdef HAVE_MPI
/// @brief A Standard timer class using the MPI MPI_Wtime
class MPITimer final : public ITimer
{
  public:
    friend class Timer;
    ~MPITimer();
  private:
    MPITimer();
    MPITimer(const MPITimer&) = delete;
    MPITimer& operator=(const MPITimer&) = delete;
    /// @brief implementation of ITimer::start
    void start() override;
    /// @brief implementation of ITimer::stop
    void stop() override;
    /// @brief implementation of ITimer::elapsedTime
    double elapsedTime() const override;
    /// @brief implementation of ITimer::summary
    std::string summary() const override;

    double itsElapsedTime;
    double itsStartTime;
    double itsStopTime;
};
#endif

#ifdef _OPENMP
// @brief A Standard timer class using the OpenMP omp_get_wtime
class OpenMPTimer final : public ITimer
{
  public:
    friend class Timer;
    ~OpenMPTimer();
  private:
    OpenMPTimer();
    OpenMPTimer(const OpenMPTimer&) = delete;
    OpenMPTimer& operator=(const OpenMPTimer &) = delete;
    /// @brief implementation of ITimer::start
    void start() override;
    /// @brief implementation of ITimer::stop
    void stop() override;
    /// @brief implementation of ITimer::elapsedTime
    double elapsedTime() const override;
    /// @brief implementation of ITimer::summary
    std::string summary() const override;

    double itsElapsedTime;
    double itsStartTime;
    double itsStopTime;
};
#endif

///@brief SectionTimer contains a map of timers
class SectionTimer final
{
  public:
    /// @brief constructor
    /// @param numTimer - number of timers to create
    SectionTimer(unsigned int numTimer);
    SectionTimer(const SectionTimer&) = delete;
    SectionTimer& operator=(const SectionTimer&) = delete;
    ~SectionTimer() {}

    /// @brief - start the timer number timerNum
    /// @param timerNum - timer number
    void start(unsigned int timerNum );

    /// @brief - stop the timer number timerNum
    /// @param timerNum - timer number
    void stop(unsigned int timerNum);
    /// @brief - log the elapsed time of each of the timer in this class
    void summary() const;
    /// @brief return the total elapsed time of all the timers in this class
    double totalElapsedTime() const;
  private:
    std::map<unsigned int, std::shared_ptr<Timer>> itsTimers;

};
}
}
#endif
