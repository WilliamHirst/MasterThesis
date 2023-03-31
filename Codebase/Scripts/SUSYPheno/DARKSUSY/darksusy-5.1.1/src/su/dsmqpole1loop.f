      real*8 function dsmqpole1loop(mqmq)
      implicit none
      real*8 mqmq,dsralph31loop
      real*8 pi
      parameter (pi=3.141592653589793238d0)
c
      dsmqpole1loop=mqmq*(1.d0+4.d0/3.d0/pi*dsralph31loop(mqmq))
      return
      end
