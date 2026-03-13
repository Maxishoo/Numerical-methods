from dataclasses import dataclass
from typing import Optional


@dataclass
class OrderRow:
    order_id: int
    created_at: str
    customer_name: str
    city_name: Optional[str]
    items_count: int
    total_amount: float
    last_payment_status: Optional[str]


class OrdersReportService:
    def __init__(
        self,
        orders_repo,
        customers_repo,
        addresses_repo,
        cities_repo,
        order_items_repo,
        payments_repo,
    ):
        self.orders_repo = orders_repo
        self.customers_repo = customers_repo
        self.addresses_repo = addresses_repo
        self.cities_repo = cities_repo
        self.order_items_repo = order_items_repo
        self.payments_repo = payments_repo


        async def build_report(self, date_from: str, date_to: str) -> list[OrderRow]:
            orders = await self.orders_repo.find_by_period(date_from, date_to)
            # SELECT * FROM orders WHERE created_at BETWEEN ? AND ?

            rows: list[OrderRow] = []

            for o in orders:
                customer = await self.customers_repo.get_by_id(o.customer_id)
            # SELECT * FROM customers WHERE id = ?

            address = await self.addresses_repo.get_by_customer_id(customer.id)
            # SELECT * FROM addresses WHERE customer_id = ? LIMIT 1

            city_name = None
            if address is not None:
                city = await self.cities_repo.get_by_id(address.city_id)
                # SELECT * FROM cities WHERE id = ?
                city_name = city.name if city else None

            items_count = await self.order_items_repo.count_by_order_id(o.id)
            # SELECT COUNT(*) FROM order_items WHERE order_id = ?

            total_amount = await self.order_items_repo.sum_total_by_order_id(o.id)
            # SELECT SUM(price * qty) FROM order_items WHERE order_id = ?

            last_payment = await self.payments_repo.get_last_by_order_id(o.id)
            # SELECT * FROM payments WHERE order_id = ? ORDER BY created_at DESC LIMIT 1

            rows.append(
                OrderRow(
                    order_id=o.id,
                    created_at=o.created_at,
                    customer_name=customer.full_name,
                    city_name=city_name,
                    items_count=items_count,
                    total_amount=float(total_amount or 0),
                    last_payment_status=last_payment.status if last_payment else None,
                )
            )

            return rows
