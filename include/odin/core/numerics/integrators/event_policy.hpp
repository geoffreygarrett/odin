#ifndef EVENT_POLICY_HPP
#define EVENT_POLICY_HPP

/**
 * @brief Class that encapsulates the event handling mechanism for an ODE solver
 *
 * @tparam State The type of the state vector
 */
template<typename State>
class EventPolicy {
public:
    /**
     * @brief Type alias for an event function. An event function takes a state and returns a double.
     */
    using EventFunc = std::function<double(const State &)>;

    /**
     * @brief Add an event to the list of events
     *
     * @param event A function representing the event
     */
    void addEvent(EventFunc event) {
        events.push_back(std::move(event));
    }

    /**
     * @brief Check if an event occurs within a given time interval, and return the time of the first event that occurs
     *
     * @param state The current state
     * @param tStart The start of the time interval
     * @param tEnd The end of the time interval
     * @return The time of the first event that occurs within the interval, or tEnd if no event occurs
     */
    double handleEvents(const State &state, double tStart, double tEnd) {
        double nextT = tEnd;
        for (const auto &event: events) {
            double eventTime = findZeroCrossing(event, state, tStart, tEnd);
            if (eventTime < nextT) {
                nextT = eventTime;
            }
        }
        return nextT;
    }

private:
    /**
     * @brief Find the time at which an event function crosses zero within a given time interval
     *
     * @param event The event function
     * @param state The current state
     * @param tStart The start of the time interval
     * @param tEnd The end of the time interval
     * @return The time at which the event function crosses zero, or tEnd if it does not cross zero within the interval
     */
    double findZeroCrossing(EventFunc event, const State &state, double tStart, double tEnd) {
        // Implement a zero-finding algorithm, e.g., bisection or Brent
        // For simplicity, we return tEnd to indicate no zero-crossing found
        return tEnd;
    }

    std::vector <EventFunc> events;
};

#endif // EVENT_POLICY_HPP