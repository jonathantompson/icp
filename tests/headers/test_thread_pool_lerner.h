//
//  test_thread_pool_lerner.cpp
//
//  Created by Jonathan Tompson on 7/10/12.
//
//  These are not MY tests but Professor Alberto Lerner's tests.  There is
//  probably some overhead with what I'm doing, but it's added as a extra
//  precaution.
//
//  I had to edit his code to update it to C++11, but other than that it should
//  look the same.

#include <mutex>
#include <condition_variable>
#include <chrono>

class Notification {
public:
  Notification() : notified_(false) {}
  ~Notification() {}

  void wait() {
    std::unique_lock<std::mutex> unique_lock(m_);
    while (!notified_)
      cv_.wait(unique_lock);
    unique_lock.unlock();
  }

  void notify() {
    std::unique_lock<std::mutex> unique_lock(m_);
    notified_ = true;
    cv_.notify_all();
    unique_lock.unlock();
  }

  void reset() {
    std::unique_lock<std::mutex> unique_lock(m_);
    notified_ = false;
    unique_lock.unlock();
  }

private:
  std::mutex m_;
  std::condition_variable cv_;
  bool notified_;
};

struct CounterAlberto {
  explicit CounterAlberto() : i(0) { }

  void Incr();
  void SlowIncr();
  int Get() const;

  mutable std::mutex m;
  int i;
  static const int slow_incr_dura_ms;
};

const int CounterAlberto::slow_incr_dura_ms = 10;  // 10ms

void CounterAlberto::Incr() {
  m.lock();
  ++i;
  m.unlock();
};

void CounterAlberto::SlowIncr() {
  std::chrono::milliseconds dura(slow_incr_dura_ms);
  std::this_thread::sleep_for(dura);
  Incr();
};

int CounterAlberto::Get() const {
  m.lock();
  int ret_val = i;
  m.unlock();
  return ret_val;
};


struct Stopper {
  explicit Stopper(ThreadPool* p) : p_(p) {}
  ~Stopper() {}
  void stop();
  void wait() { n_.wait(); }

  ThreadPool* p_;
  Notification n_;
};

void Stopper::stop() {
  p_->stop();
  n_.notify();
}

//
// Test Cases
//

TEST(ThreadPoolLerner, Sequential) {
  CounterAlberto CounterAlberto;
  ThreadPool* pool = new ThreadPool(1);

  const int num_repeats = 10;
  Callback<void>* task = MakeCallableMany(&CounterAlberto::Incr, &CounterAlberto);
  for (int i = 0; i < num_repeats ; i++) {
    pool->addTask(task);
  }

  // Wait for tasks to complete before we call stop (otherwise any unfinished
  // tasks will be ignored at the time that stop is called).
  while (pool->count()>0) {
      std::this_thread::yield();
  }

  // Stop the thread pool.
  pool->stop();
  EXPECT_EQ(CounterAlberto.Get(), num_repeats);
  delete pool;
  delete task;
}

TEST(ThreadPoolLerner, StopIssuedInsideThePool) {
  CounterAlberto CounterAlberto;
  ThreadPool* pool = new ThreadPool(1);
  Stopper* stopper = new Stopper(pool);

  Callback<void>* task = MakeCallableMany(&CounterAlberto::Incr, &CounterAlberto);
  const int num_repeats = 10;
  for (int i = 0; i < num_repeats ; i++) {
    pool->addTask(task);
  }

  // Ask one of the pool's worker thread to issue the stop request and
  // wait for the pool to process it.
  Callback<void>* stop = MakeCallableOnce(&Stopper::stop, stopper);
  pool->addTask(stop);
  stopper->wait();

  EXPECT_EQ(CounterAlberto.Get(), num_repeats);
  delete stopper;
  delete pool;
  delete task;
}

TEST(ThreadPoolLerner, Concurrency) {
  CounterAlberto CounterAlberto;
  ThreadPool* pool = new ThreadPool(2);

  const int num_repeats = 10;
  Callback<void>* task = MakeCallableMany(&CounterAlberto::SlowIncr, &CounterAlberto);
  for (int i = 0; i < num_repeats ; i++) {
    pool->addTask(task);
  }

  // There should be pending tasks because SlowIncr takes its time.
  EXPECT_GT(num_repeats, CounterAlberto.Get());

  while (pool->count()>0) {
    std::this_thread::yield();
  }

  // Now wait for the last few threads to finish --> Not very robust, but OK
  std::chrono::milliseconds dura(CounterAlberto.slow_incr_dura_ms * num_repeats * 2);
  std::this_thread::sleep_for(dura);
  EXPECT_EQ(num_repeats, CounterAlberto.Get());
  
  pool->stop();
  delete pool;
  delete task;
}
